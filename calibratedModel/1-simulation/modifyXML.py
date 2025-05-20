import os
from lxml import etree as ET

file = "./calibrationSimulation.xml"
output = "./modified.xml"

context = ET.iterparse(file, events=("start", "end"))
root = None

for event, elem in context:
    if event == "start" and elem.tag == "alignment":
        elem.clear()  # remove the alignment
    elif root is None:
        root = elem  
# Step 1: add samplefromprior
# for mcmc in root.findall(".//run"):
#     if mcmc.get("id") == "MCMC" and mcmc.get("spec") == "MCMC":
#         mcmc.set("sampleFromPrior", "true")
# Step 2: delete unwanted operators
for operator in root.findall(".//operator"):
    operator_id = operator.get("id", "")
    if operator_id.startswith("deltaExchangePair") or operator_id.startswith("swapPair"):
        parent = operator.getparent()
        if parent is not None:
            parent.remove(operator)

# Step 3: add new operator
first_operator = root.find(".//operator")
if first_operator is not None:
    parent = first_operator.getparent()
    if parent is not None:
        
        if parent.find(".//operator[@id='deltaExchangeHomozygous']") is None:
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

# Step 4: change log file name
for logFile in root.findall(".//logger[@id='Logger']"):
    logFile.set("fileName", "treesFromPrior.log")

for treeFile in root.findall(".//logger[@id='tree.treeLogger']"):
    treeFile.set("fileName", "treesFromPrior.trees")

# Step 5: save the modified xml
tree = ET.ElementTree(root)
tree.write(output, pretty_print=True, xml_declaration=True, encoding="UTF-8")

print(f"Modified XML saved as {output}")
