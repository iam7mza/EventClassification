#!/bin/bash
# filepath: /home/hamza/ML/script.sh

# Set the sample size as a variable
SAMPLE_SIZE=1000

# Create output directory
mkdir -p csv_output

# Loop over all .root files in the current directory
for file in *.root; do
    if [ -f "$file" ]; then
        echo "Processing: $file with sample size $SAMPLE_SIZE"
        /home/hamza/ML/rootToCSV "$file" $SAMPLE_SIZE
        
        # Move the generated CSV to csv_output directory
        csv_file="${file%.root}_sample${SAMPLE_SIZE}.csv"
        if [ -f "$csv_file" ]; then
            mv "$csv_file" csv_output/
            echo "  -> Moved to csv_output/$csv_file"
        fi
    fi
done

echo "Done! CSV files are in csv_output/"
