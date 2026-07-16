#!/bin/bash

# Script to copy files from ../3-xml and modify lphy files
# Preserves directory structure

set -e  # Exit on error

SOURCE_DIR="../3-xml"
CURRENT_DIR="."

echo "Starting copy and modification process..."
echo "Source directory: $SOURCE_DIR"
echo "Target directory: $CURRENT_DIR"
echo ""

# Step 1: Copy all .lphy files from ../3-xml/*Scripts/ while replacing 2-variantCalling with 2-variantCalling-method2
echo "Step 1: Copying and modifying .lphy files from $SOURCE_DIR/*Scripts/..."
find "$SOURCE_DIR"/*Scripts -name ".ipynb_checkpoints" -prune -o -name "*.lphy" -type f -print 2>/dev/null | while read lphy_file; do
    # Create relative path
    relative_path="${lphy_file#$SOURCE_DIR/}"
    target_file="$CURRENT_DIR/$relative_path"
    target_dir=$(dirname "$target_file")
    
    # Create directory if it doesn't exist
    mkdir -p "$target_dir"
    
    # Copy file and replace 2-variantCalling with 2-variantCalling-method2
    sed 's/2-variantCalling\b/2-variantCalling-method2/g' "$lphy_file" > "$target_file"
    echo "Copied and modified: $relative_path"
done

# Step 2: Copy all *-* files from ../3-xml while preserving directory structure
echo ""
echo "Step 2: Copying *-* files from $SOURCE_DIR..."
find "$SOURCE_DIR" -name ".ipynb_checkpoints" -prune -o -name "*-*" -type f -print 2>/dev/null | while read source_file; do
    # Skip directories and only copy regular files
    if [ -f "$source_file" ]; then
        # Create relative path
        relative_path="${source_file#$SOURCE_DIR/}"
        target_file="$CURRENT_DIR/$relative_path"
        target_dir=$(dirname "$target_file")
        
        # Create directory if it doesn't exist
        mkdir -p "$target_dir"
        
        # Copy file (skip if it already exists from lphy copy)
        if [ ! -f "$target_file" ]; then
            cp "$source_file" "$target_file"
            echo "Copied: $relative_path"
        else
            echo "Skipped (already exists): $relative_path"
        fi
    fi
done

echo ""
echo "✓ Process completed successfully!"
echo ""
echo "Summary:"
echo "- Copied and modified .lphy files from $SOURCE_DIR/*Scripts/"
echo "  (2-variantCalling replaced with 2-variantCalling-method2)"
echo "- Copied *-* files from $SOURCE_DIR"
echo "- Preserved directory structure throughout"
echo "- Excluded .ipynb_checkpoints folders"