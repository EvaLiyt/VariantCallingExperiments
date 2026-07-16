import os
import shutil
from lxml import etree as ET

# Define a base output directory
output_folder = "./xmls"
os.makedirs(output_folder, exist_ok=True)

# Clear the output folder
# for filename in os.listdir(output_folder):
#     file_path = os.path.join(output_folder, filename)
#     if os.path.isfile(file_path):
#         os.remove(file_path)
#     elif os.path.isdir(file_path):
#         shutil.rmtree(file_path)

def modify_xml_files(input_folder, output_folder):
    # Process each XML file in the input folder
    for filename in os.listdir(input_folder):
        if filename.endswith(".xml"):
            filepath = os.path.join(input_folder, filename)
            xml = ET.parse(filepath)
            root = xml.getroot()

            # Step 2: Remove specific <operator> nodes
            for operator in root.findall(".//operator"):
                operator_id = operator.get("id", "")
                if operator_id.startswith("deltaExchangePair") or operator_id.startswith("swapPair"):
                    parent = operator.getparent()
                    if parent is not None:
                        parent.remove(operator)

            # Step 3: Add new operators before the first <operator>
            first_operator = root.find(".//operator")
            if first_operator is not None:
                parent = first_operator.getparent()
                if parent is not None:
                    # Avoid duplication
                    if parent.find(".//operator[@id='deltaExchangeHomozygous']") is None:
                        # Homozygous
                        homo_operator = ET.Element(
                            "operator",
                            attrib={
                                "id": "deltaExchangeHomozygous",
                                "spec": "operator.DeltaExchangeOperator",
                                "delta": "0.0625",
                                "weight": "4.0",
                            }
                        )
                        ET.SubElement(homo_operator, "parameter", attrib={"idref": "pi"})
                        ET.SubElement(homo_operator, "weightvector", attrib={
                            "id": "IntegerParameter",
                            "spec": "parameter.IntegerParameter",
                            "dimension": "16",
                            "lower": "0",
                            "upper": "0"
                        }).text = "1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1"
                        parent.insert(parent.index(first_operator), homo_operator)

                        # Heterozygous
                        hetero_operator = ET.Element(
                            "operator",
                            attrib={
                                "id": "deltaExchangeHeterozygous",
                                "spec": "operator.DeltaExchangeOperator",
                                "delta": "0.0625",
                                "weight": "12.0",
                            }
                        )
                        ET.SubElement(hetero_operator, "parameter", attrib={"idref": "pi"})
                        ET.SubElement(hetero_operator, "weightvector", attrib={
                            "id": "IntegerParameter1",
                            "spec": "parameter.IntegerParameter",
                            "dimension": "16",
                            "lower": "0",
                            "upper": "0"
                        }).text = "0 1 1 1 1 0 1 1 1 1 0 1 1 1 1 0"
                        parent.insert(parent.index(first_operator), hetero_operator)
            else:
                print(f"Warning: No <operator> elements found in: {filename}")

            output_filename = f"{filename}"
            output_path = os.path.join(output_folder, output_filename)
            xml.write(output_path, pretty_print=True, xml_declaration=True, encoding="UTF-8")
            print(f"Processed and saved: {output_path}")

# List of folders to process
directories = [
    "./ErrModelScripts",
    "./wtErrModelScripts",
    "./TrueErrModelScripts",
    "./TruewtErrModelScripts"
]

# Modify all XML files
for folder in directories:
    modify_xml_files(input_folder=folder, output_folder=output_folder)
