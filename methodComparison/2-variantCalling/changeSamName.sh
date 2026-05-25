#!/bin/bash

# Fix SAM files:
#   1. Replace all sequence names {number}-{1,2} with chr1 everywhere
#   2. Keep only the FIRST @SQ line, remove all subsequent @SQ lines
#
# Usage: bash fix_sam_headers.sh
# Run from the PARENT folder of bulkSeqResult/
# Edits .sam files IN-PLACE — back up first if needed!

SAM_DIR="bulkSeqResult"

if [ ! -d "$SAM_DIR" ]; then
  echo "ERROR: Directory '$SAM_DIR' not found. Run this from its parent folder."
  exit 1
fi

SAM_FILES=("$SAM_DIR"/*.sam)

if [ ! -e "${SAM_FILES[0]}" ]; then
  echo "No .sam files found in '$SAM_DIR'."
  exit 1
fi

echo "Found ${#SAM_FILES[@]} SAM file(s). Processing..."

for sam in "${SAM_FILES[@]}"; do
  # Step 1: Replace all {number}-{1,2} sequence names with chr1
  sed -i -E 's/chr1-[12]/chr1/g' "$sam"

  # Step 2: Remove duplicate @SQ lines — keep only the first occurrence
  awk '/^@SQ/ { if (seen++) next } 1' "$sam" > "${sam}.tmp" && mv "${sam}.tmp" "$sam"

  echo "  Updated: $sam"
done

echo "Done."
