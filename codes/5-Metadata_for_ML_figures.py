import csv
from Bio import Phylo

print("--- Generating Metadata Dictionary for Annotation ---")

tree_file = "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/lep_phylogeny_final.treefile"
csv_out = "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/leprosy_metadata.csv"

# Load the tree
tree = Phylo.read(tree_file, "newick")

# Create the CSV file
with open(csv_out, mode='w', newline='') as file:
    writer = csv.writer(file)
    # Write the header row
    writer.writerow(["ID", "Host", "Era", "Location", "Color_Code"])

    # Loop through every sample in the tree
    for leaf in tree.get_terminals():
        # Clean the ID to match the exact names in your visualizer
        node_id = leaf.name.split('/')[-1].replace('.final.bam', '')

        # Skip the M. lepromatosis outgroups so they don't clutter your data
        if tree.distance(leaf) > 0.05:
            continue

        # 1. The Animal Reservoir (ENA Data)
        if "SRR367" in node_id:
            writer.writerow([node_id, "Red Squirrel", "Modern", "United Kingdom", "#FF0000"]) # Red

        # 2. The Historical Epidemic (PRJNA200950)
        elif "SRR847" in node_id:
            writer.writerow([node_id, "Human", "Medieval (10th-12th C.)", "Winchester, UK", "#0000FF"]) # Blue

        # 3. The Global Outgroups (PRJNA317287 & others)
        else:
            writer.writerow([node_id, "Human", "Modern", "Global", "#808080"]) # Gray

print(f"✅ Metadata dictionary successfully saved to: {csv_out}")
