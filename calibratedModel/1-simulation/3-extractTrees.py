import json
import random
import shutil
import os
import re

# Load configuration
with open("extractConfig.json", "r") as f:
    config = json.load(f)

output_path = "./trueTrees"

# Remove old output directory
if os.path.exists(output_path):
    shutil.rmtree(output_path)
    print(f"The directory {output_path} and all its contents were deleted.")

os.makedirs(output_path)
print(f"The directory {output_path} was created.")

tree_file = config["tree_file"]
num_trees = config["num_trees"]
output_prefix = config["output_prefix"]
output_suffix = config["output_suffix"]

# Read the whole file
with open(tree_file, "r") as f:
    lines = f.readlines()

# Parse taxa block and Translate mapping
taxa_block = []
translate_dict = {}
in_taxa = False
in_translate = False
header_end_index = 0

for i, line in enumerate(lines):
    stripped = line.strip()
    if stripped.lower().startswith("begin taxa"):
        in_taxa = True
    elif stripped.lower().startswith("end;") and in_taxa:
        in_taxa = False
        taxa_block.append("End;\n")
        continue
    elif in_taxa:
        taxa_block.append(line)
    
    if stripped.lower().startswith("translate"):
        in_translate = True
    elif in_translate:
        if ";" in stripped:
            in_translate = False
            stripped = stripped.replace(";", "")
        stripped = stripped.strip().rstrip(",")
        if stripped:
            parts = stripped.split()
            if len(parts) == 2:
                index, name = parts
                translate_dict[index] = name

    if stripped.startswith("tree "):
        header_end_index = i
        break

# Extract all tree lines
trees = [line.strip() for line in lines if line.strip().startswith("tree ")]

# Sample trees
selected_trees = random.sample(trees, num_trees)

# Function to replace only tip numbers before colons
def replace_tree_taxa(tree_line, translate_dict):
    # Match only numbers that come *before* a colon (not affecting branch lengths)
    def repl(match):
        num = match.group(1)
        return translate_dict.get(num, num) + ":"
    
    return re.sub(r'\b(\d+):', repl, tree_line)

# Output selected trees with restored taxa names
for i, tree in enumerate(selected_trees):
    output_file = f"{output_path}/{output_prefix}{i}{output_suffix}"
    with open(output_file, "w") as f:
        # Write only NEXUS + taxa block (skip translate)
        f.write("#NEXUS\n\n")
        f.writelines(taxa_block)
        f.write("\nBegin trees;\n")
        restored_tree = replace_tree_taxa(tree, translate_dict)
        f.write(restored_tree + "\n")
        f.write("End;\n")
