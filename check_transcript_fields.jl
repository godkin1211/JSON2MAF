using JSON2MAF
using JSON3
using CodecZlib

println("Loading JSON file...")
json_data = open("P.hard-filtered.vcf.annotated.json.gz", "r") do file
    stream = GzipDecompressorStream(file)
    JSON3.read(stream)
end

println("Analyzing transcript structure...\n")

if haskey(json_data, :positions) && length(json_data.positions) > 0
    # 檢查前幾個 position
    for pos_idx in 1:min(3, length(json_data.positions))
        pos = json_data.positions[pos_idx]

        if haskey(pos, :transcripts) && length(pos.transcripts) > 0
            println("="^60)
            println("Position $pos_idx: $(pos.chromosome):$(pos.start)")
            println("="^60)

            # 檢查前幾個 transcript
            for (i, trans) in enumerate(pos.transcripts[1:min(3, length(pos.transcripts))])
                println("\nTranscript $i:")
                println("  ID: ", get(trans, :transcript, "N/A"))
                println("  Gene: ", get(trans, :geneSymbol, get(trans, :gene, "N/A")))
                println("  Source: ", get(trans, :source, "N/A"))
                println("  Biotype: ", get(trans, :bioType, "N/A"))

                # 檢查可能的 MANE/canonical 欄位
                mane_fields = [:isCanonical, :isManeSelect, :maneSelect, :canonical,
                              :MANE, :mane, :isManePlus, :manePlus]

                found_mane = false
                for field in mane_fields
                    if haskey(trans, field)
                        println("  ✓ $field = ", get(trans, field, nothing))
                        found_mane = true
                    end
                end

                if !found_mane
                    println("  (No MANE/canonical markers found)")
                end
            end
            println()
        end
    end

    # 統計有多少 transcript 包含 MANE/canonical 標記
    println("\n" * "="^60)
    println("Statistics across all positions:")
    println("="^60)

    total_transcripts = 0
    has_canonical = 0
    has_mane = 0

    for pos in json_data.positions
        if haskey(pos, :transcripts)
            for trans in pos.transcripts
                total_transcripts += 1
                if haskey(trans, :isCanonical) || haskey(trans, :canonical)
                    has_canonical += 1
                end
                if haskey(trans, :isManeSelect) || haskey(trans, :maneSelect) ||
                   haskey(trans, :MANE) || haskey(trans, :mane)
                    has_mane += 1
                end
            end
        end
    end

    println("Total transcripts: $total_transcripts")
    println("With canonical marker: $has_canonical ($(round(100*has_canonical/total_transcripts, digits=2))%)")
    println("With MANE marker: $has_mane ($(round(100*has_mane/total_transcripts, digits=2))%)")

    # 列出所有unique的transcript欄位
    println("\n" * "="^60)
    println("All unique transcript field names found:")
    println("="^60)

    all_fields = Set{Symbol}()
    for pos in json_data.positions[1:min(100, length(json_data.positions))]
        if haskey(pos, :transcripts)
            for trans in pos.transcripts
                union!(all_fields, keys(trans))
            end
        end
    end

    for field in sort(collect(all_fields))
        println("  - $field")
    end
else
    println("No positions found in JSON")
end
