#!/bin/bash

# Usage: ./find_large_dirs_and_files.sh /source/path 100M

SOURCE_DIR="$1"
SIZE_LIMIT="$2"

if [ -z "$SOURCE_DIR" ] || [ -z "$SIZE_LIMIT" ]; then
    echo "Usage: $0 <source_dir> <size_limit>"
    echo "Example: $0 /home/user/data 500M"
    exit 1
fi

# Convert human-readable size to bytes
SIZE_LIMIT_BYTES=$(numfmt --from=iec "$SIZE_LIMIT")

echo "Scanning directories in: $SOURCE_DIR"
echo "Size threshold: $SIZE_LIMIT ($SIZE_LIMIT_BYTES bytes)"
echo "-------------------------------------------------------"

# Check each top-level directory
find "$SOURCE_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r dir; do
    dir_size=$(du -s --block-size=1 "$dir" | cut -f1)

    if [ "$dir_size" -ge "$SIZE_LIMIT_BYTES" ]; then
        echo ""
        echo "📁 Directory over threshold: $dir (Size: $(numfmt --to=iec $dir_size))"

        # Now find large files inside this directory
        echo "   🔍 Large files inside:"
        find "$dir" -type f -size +"$SIZE_LIMIT" -exec du -h {} + | sort -hr
    fi
done