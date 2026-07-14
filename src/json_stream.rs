/// Single-pass streaming reader for Nirvana-style JSON documents of the form
/// `{ "header": <H>, "positions": [ <T>, <T>, ... ], ... }`.
///
/// The whole point of this module: never materialize the document (or even
/// the whole `positions` array) as a `serde_json::Value` DOM or a `Vec<T>`.
/// Each element of `positions` is deserialized and hand off to a caller
/// callback one at a time, then dropped, so peak memory stays proportional to
/// a single position rather than the whole (potentially multi-GB decompressed)
/// file. This is what makes JSON2MAF able to process large DRAGEN/Nirvana
/// annotation files on machines with limited RAM.
use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use serde::de::{DeserializeOwned, DeserializeSeed, Error as DeError, IgnoredAny, MapAccess, SeqAccess, Visitor};
use serde::Deserializer as _;
use std::fmt;
use std::fs::File;
use std::io::BufReader;

/// Streams `positions` array elements of type `T` out of the gzipped JSON
/// document at `file_path`, calling `on_header` once (as soon as the
/// `header` key is parsed) and `on_item` once per position element, in file
/// order. Returns the parsed header.
///
/// Nirvana always emits `header` before `positions`; `on_header` is
/// guaranteed to have already run by the time `on_item` is first called.
pub fn stream_positions<T, H, OnHeader, OnItem>(
    file_path: &str,
    on_header: OnHeader,
    on_item: OnItem,
) -> Result<H>
where
    T: DeserializeOwned,
    H: DeserializeOwned,
    OnHeader: FnMut(&H) -> Result<()>,
    OnItem: FnMut(T) -> Result<()>,
{
    let file = File::open(file_path)
        .with_context(|| format!("Failed to open input file: {}", file_path))?;
    let reader = BufReader::new(file);
    let decoder = MultiGzDecoder::new(reader);
    let buffered = BufReader::with_capacity(1024 * 1024, decoder);
    let mut json_de = serde_json::Deserializer::from_reader(buffered);

    let visitor = DocumentVisitor {
        on_header,
        on_item,
        _marker: std::marker::PhantomData::<T>,
        _h: std::marker::PhantomData::<H>,
    };

    let header = json_de
        .deserialize_map(visitor)
        .context("Failed to parse Nirvana JSON document")?;

    json_de
        .end()
        .context("Unexpected trailing data after JSON document")?;

    Ok(header)
}

struct DocumentVisitor<T, H, OnHeader, OnItem> {
    on_header: OnHeader,
    on_item: OnItem,
    _marker: std::marker::PhantomData<T>,
    #[allow(dead_code)]
    _h: std::marker::PhantomData<H>,
}

impl<'de, T, H, OnHeader, OnItem> Visitor<'de> for DocumentVisitor<T, H, OnHeader, OnItem>
where
    T: DeserializeOwned,
    H: DeserializeOwned,
    OnHeader: FnMut(&H) -> Result<()>,
    OnItem: FnMut(T) -> Result<()>,
{
    type Value = H;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("a Nirvana JSON document object with \"header\" and \"positions\"")
    }

    fn visit_map<A>(mut self, mut map: A) -> std::result::Result<Self::Value, A::Error>
    where
        A: MapAccess<'de>,
    {
        let mut header: Option<H> = None;

        while let Some(key) = map.next_key::<String>()? {
            match key.as_str() {
                "header" => {
                    let parsed: H = map.next_value()?;
                    (self.on_header)(&parsed).map_err(DeError::custom)?;
                    header = Some(parsed);
                }
                "positions" => {
                    map.next_value_seed(PositionsSeed {
                        on_item: &mut self.on_item,
                        _marker: std::marker::PhantomData::<T>,
                    })?;
                }
                _ => {
                    map.next_value::<IgnoredAny>()?;
                }
            }
        }

        header.ok_or_else(|| DeError::custom("No \"header\" field found in JSON document"))
    }
}

struct PositionsSeed<'f, T, OnItem> {
    on_item: &'f mut OnItem,
    _marker: std::marker::PhantomData<T>,
}

impl<'de, 'f, T, OnItem> DeserializeSeed<'de> for PositionsSeed<'f, T, OnItem>
where
    T: DeserializeOwned,
    OnItem: FnMut(T) -> Result<()>,
{
    type Value = ();

    fn deserialize<D>(self, deserializer: D) -> std::result::Result<Self::Value, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        deserializer.deserialize_seq(self)
    }
}

impl<'de, 'f, T, OnItem> Visitor<'de> for PositionsSeed<'f, T, OnItem>
where
    T: DeserializeOwned,
    OnItem: FnMut(T) -> Result<()>,
{
    type Value = ();

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("an array of position objects")
    }

    fn visit_seq<A>(self, mut seq: A) -> std::result::Result<Self::Value, A::Error>
    where
        A: SeqAccess<'de>,
    {
        while let Some(item) = seq.next_element::<T>()? {
            (self.on_item)(item).map_err(DeError::custom)?;
        }
        Ok(())
    }
}
