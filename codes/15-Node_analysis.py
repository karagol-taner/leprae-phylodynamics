from Bio import Phylo
import re

print("--- Advanced Bayesian Node Analysis ---")
tree_file = "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/lep_MCC_dated.tree"
tree = next(Phylo.parse(tree_file, "nexus"))

# 1. Isolate the UK Squirrel Clade explicitly
# (Grabbing known UK SRRs to bypass the ancient lepromatosis outgroup)
target_sq_ids = ["SRR3672759", "SRR3672758", "SRR3672757"]
uk_sq_tips = [t for t in tree.get_terminals() if any(sq_id in t.name for sq_id in target_sq_ids)]
sq_mrca = tree.common_ancestor(uk_sq_tips)

# 2. Get the Parent Node (The Divergence from the Winchester Sister Lineage)
path_to_sq = tree.get_path(sq_mrca)
divergence_node = path_to_sq[-2] if len(path_to_sq) > 1 else tree.root

def extract_dates(node):
    comment = node.comment
    if not comment: return None
    med = re.search(r'height_median=([\d\.E-]+)', comment)
    hpd = re.search(r'height_95%_HPD=\{([\d\.E-]+),([\d\.E-]+)\}', comment)
    post = re.search(r'posterior=([\d\.E-]+)', comment)
    if med and hpd:
        return {
            "median": 2020 - float(med.group(1)),
            "lower": 2020 - float(hpd.group(2)), # HPD upper height = oldest year
            "upper": 2020 - float(hpd.group(1)), # HPD lower height = youngest year
            "posterior": float(post.group(1)) if post else 1.0
        }
    return None

sq_dates = extract_dates(sq_mrca)
div_dates = extract_dates(divergence_node)

if sq_dates and div_dates:
    print(f"1. Sister-Branch Divergence Node (When Squirrels split from Humans):")
    print(f"   Median: {div_dates['median']:.0f} AD (95% HPD: {div_dates['lower']:.0f} - {div_dates['upper']:.0f})")
    print(f"   Posterior Support: {div_dates['posterior']}")

    print(f"\n2. Squirrel Clade MRCA (When the squirrel reservoir began radiating):")
    print(f"   Median: {sq_dates['median']:.0f} AD (95% HPD: {sq_dates['lower']:.0f} - {sq_dates['upper']:.0f})")
    print(f"   Posterior Support: {sq_dates['posterior']}")

    # Mathematical test: Evolutionary Stem Length
    stem_length = sq_dates['median'] - div_dates['median']
    print(f"\n3. Evolutionary Stem Branch Length:")
    print(f"   {stem_length:.0f} years of independent genetic drift before the squirrel clade expanded.")
else:
    print("Could not extract BEAST metadata.")
