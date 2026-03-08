#!/bin/bash

# Output file name
OUTPUT="merged_output.txt"

# Clear output file if it exists
> "$OUTPUT"

echo "Merging files into $OUTPUT ..."
echo "Generated on $(date)" >> "$OUTPUT"
echo "=====================================" >> "$OUTPUT"
echo "" >> "$OUTPUT"

# Loop through include and src folders
for dir in include src; do
    if [ -d "$dir" ]; then
        for file in "$dir"/*.{cuh,hpp,cpp,cu}; do
            if [ -f "$file" ]; then
                echo "===== FILE: $file =====" >> "$OUTPUT"
                cat "$file" >> "$OUTPUT"
                echo "" >> "$OUTPUT"
                echo "" >> "$OUTPUT"
            fi
        done
    fi
done

echo "Done. Output written to $OUTPUT"
