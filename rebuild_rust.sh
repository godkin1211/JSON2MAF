#!/bin/bash
# Force complete rebuild of Rust version

set -e

echo "=========================================="
echo "Force Rebuilding Rust JSON2MAF"
echo "=========================================="
echo ""

cd "$(dirname "$0")/rust_version"

echo "Step 1: Removing old binary..."
rm -f target/release/json2maf
echo "✓ Old binary removed"
echo ""

echo "Step 2: Cleaning build cache..."
cargo clean
echo "✓ Build cache cleaned"
echo ""

echo "Step 3: Removing lock file..."
rm -f Cargo.lock
echo "✓ Lock file removed"
echo ""

echo "Step 4: Building release binary..."
cargo build --release
echo "✓ Build complete"
echo ""

echo "Step 5: Verifying binary..."
if [ -f "target/release/json2maf" ]; then
    echo "✓ Binary exists at: $(pwd)/target/release/json2maf"
    echo "  Size: $(ls -lh target/release/json2maf | awk '{print $5}')"
    echo "  Timestamp: $(stat -c '%y' target/release/json2maf)"
else
    echo "✗ ERROR: Binary not found!"
    exit 1
fi
echo ""

echo "Step 6: Verifying code has MultiGzDecoder..."
if grep -q "MultiGzDecoder" src/parser.rs; then
    echo "✓ Code uses MultiGzDecoder (correct)"
else
    echo "✗ ERROR: Code does not use MultiGzDecoder!"
    exit 1
fi
echo ""

echo "=========================================="
echo "Rebuild Complete!"
echo "=========================================="
echo ""
echo "To test, run:"
echo "  ./rust_version/target/release/json2maf -i P.hard-filtered.vcf.annotated.json.gz -o output.maf --verbose"
