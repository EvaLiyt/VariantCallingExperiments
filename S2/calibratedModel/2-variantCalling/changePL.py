import os

input_folder = "./merge"
output_folder = "./merge"

# Loop through all files in the input folder
for filename in os.listdir(input_folder):
    # Process files starting with "merged_" and ending with ".vcf"
    if filename.startswith("merged_") and filename.endswith(".vcf"):
        input_path = os.path.join(input_folder, filename)
        
        # Extract the r_value (e.g., r0, r1) from the filename
        r_value = filename.replace("merged_", "").replace(".vcf", "")

        # Construct the output filename as variants_r0.vcf
        output_filename = f"variants_{r_value}.vcf"
        output_path = os.path.join(output_folder, output_filename)

        # Check if the output file already exists, and delete it if it does
        if os.path.exists(output_path):
            print(f"{output_path} exists. Overwriting...")
            os.remove(output_path)

        # Open the input file and output file for reading and writing respectively
        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                # If the line is the FORMAT header with ID=PL, modify it
                if line.startswith("##FORMAT=<ID=PL"):
                    line = line.replace("Number=1", "Number=G")

                # If the line starts with "#CHROM", process cell names
                if line.startswith("#CHROM"):
                    parts = line.strip().split('\t')
                    header = parts[:9]
                    cells = [cell.split('_', 1)[1] if '_' in cell else cell for cell in parts[9:]]
                    line = '\t'.join(header + cells) + '\n'


                # Write the (modified or unmodified) line to the output file
                outfile.write(line)

        print(f"Processed {input_path} and saved to {output_path}")

print("Finished modifying VCF files.")
