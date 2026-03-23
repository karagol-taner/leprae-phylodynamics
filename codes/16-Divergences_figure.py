import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from collections import defaultdict

# ── CONFIG ────────────────────────────────────────────────────
# ── CONFIG ────────────────────────────────────────────────────
TREE_FILE = "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/lep_MCC_dated.tree"
META_FILE = "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/leprosy_metadata.csv"
OUT_FIG   = "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/Figure_Divergence_Analysis_Final_c.png"

REFERENCE_YEAR = 2020
START_YEAR = 400
FUR_TRADE_START, FUR_TRADE_END = 900, 1300

# ── 1. TREE PARSER ─────────────────────────────────────
def parse_nexus_tree(tree_str):
    tokens = re.findall(r'\(|\)|,|\[&.*?\]|:[\d.E-]+|[\d\w/._-]+', tree_str)
    nodes, stack, nid_counter = {}, [], 0
    current_node = None

    for token in tokens:
        if token == '(':
            nid_counter += 1
            new_node = {'id': nid_counter, 'children': [], 'height': None, 'is_tip': False, 'label': None}
            if stack: stack[-1]['children'].append(new_node)
            stack.append(new_node)
        elif token == ')':
            current_node = stack.pop()
            nodes[current_node['id']] = current_node
        elif token == ',': continue
        elif token.startswith('[&'):
            h_match = re.search(r'height=([0-9\.E-]+)', token)
            if h_m := h_match:
                if current_node: current_node['height'] = float(h_m.group(1))
        elif token.startswith(':'): continue
        else:
            nid_counter += 1
            tip_node = {'id': nid_counter, 'children': [], 'height': None, 'is_tip': True, 'label': token}
            if stack: stack[-1]['children'].append(tip_node)
            nodes[nid_counter] = tip_node
            current_node = tip_node
    return nodes

def get_descendant_tips(node):
    if node['is_tip']: return {node['label']}
    tips = set()
    for child in node['children']:
        tips.update(get_descendant_tips(child))
    return tips

# ── 2. DATA LOADING & MAPPING ─────────────────────────────────
with open(TREE_FILE, 'r') as f:
    content = f.read()

tree_match = re.search(r'tree\s+[^=]+=\s+(.*?;)', content, re.DOTALL | re.IGNORECASE)
nodes_map = parse_nexus_tree(tree_match.group(1))

meta = pd.read_csv(META_FILE)
id_col = meta.columns[0]
# Standardize IDs to match SRR names [cite: 5]
meta["CleanID"] = meta[id_col].astype(str).str.strip()
meta_dict = meta.set_index('CleanID').to_dict('index')

translate_dict = {}
translate_block = re.search(r'Translate\s+(.*?);', content, re.DOTALL | re.IGNORECASE)
if translate_block:
    for line in translate_block.group(1).strip().split('\n'):
        parts = line.strip().split()
        if len(parts) >= 2:
            idx, label = parts[0], parts[1].strip(',').strip("'")
            # Extract SRR ID from the file path in the tree [cite: 7, 8]
            srr_match = re.search(r'(SRR\d+)', label)
            translate_dict[idx] = srr_match.group(1) if srr_match else label

tip_info = {idx: meta_dict.get(srr, {'Host': 'Unknown', 'Era': 'Unknown'})
            for idx, srr in translate_dict.items()}

# Updated matching for "Medieval (10th-12th C.)"
modern_human_tips = {idx for idx, info in tip_info.items() if info['Era'] == 'Modern' and info['Host'] == 'Human'}
medieval_human_tips = {idx for idx, info in tip_info.items() if 'Medieval' in str(info['Era']) and info['Host'] == 'Human'}
squirrel_tips = {idx for idx, info in tip_info.items() if info['Host'] == 'Red Squirrel'}

# ── 3. EXTRACT DIVERGENCE KINETICS ───────────────────────────
modern_h_years, medieval_h_years, squirrel_years = [], [], []

for nid, node in nodes_map.items():
    if node['is_tip'] or node['height'] is None: continue
    desc_tips = get_descendant_tips(node)
    total = len(desc_tips)

    m_h_count = len(desc_tips.intersection(modern_human_tips))
    med_h_count = len(desc_tips.intersection(medieval_human_tips))
    sq_count = len(desc_tips.intersection(squirrel_tips))

    year = REFERENCE_YEAR - node['height']
    if START_YEAR <= year <= REFERENCE_YEAR:
        if sq_count > 1 and sq_count / total > 0.8:
            squirrel_years.append(year)
        elif m_h_count > 1 and m_h_count / total > 0.8:
            modern_h_years.append(year)
        # Medieval lineage detection
        elif med_h_count > 1 and med_h_count / total > 0.8:
            medieval_h_years.append(year)

# ── 4. VISUALIZATION ──────────────────────────────────────────
plt.rcParams['font.family'] = 'sans-serif'
fig, axes = plt.subplots(3, 1, figsize=(14, 15))
fig.patch.set_facecolor("white")
bins = list(range(START_YEAR, 2025, 50))

def style_panel(ax, years, color, label, title, show_fur_trade=False, post_1100_color=None):
    ax.set_facecolor("#FAFAF8")
    if show_fur_trade:
        ax.axvspan(FUR_TRADE_START, FUR_TRADE_END, alpha=0.15, color="#8B4513", label="Fur Trade Peak (900–1300 AD)")
    ax.axvline(1050, color="#8b0000", lw=1.5, ls="--", alpha=0.8, label="Winchester Epidemic ~1050 AD")

    if years:
        counts, edges = np.histogram(years, bins=bins)
        mids = [(edges[j] + edges[j+1]) / 2 for j in range(len(edges)-1)]

        # Panel C specific color shift logic
        bar_colors = [post_1100_color if (post_1100_color and m > 1100) else color for m in mids]
        ax.bar(mids, counts, width=45, color=bar_colors, alpha=0.75, edgecolor="white", label=label, zorder=3)

    ax.set_xlim(400, 2020)
    ax.set_ylim(0, max(3, (max(np.histogram(years, bins=bins)[0]) + 1) if years else 3))
    ax.tick_params(axis='both', which='major', bottom=True, left=True, labelbottom=True, direction='out', length=6)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontweight('bold')
    ax.set_ylabel("Lineage Splits", fontsize=10, fontweight="bold")
    ax.set_title(title, fontsize=13, fontweight="bold", loc="left")
    ax.set_xlabel("Calendar Year (AD)", fontsize=12, fontweight="bold", labelpad=10)
    ax.legend(fontsize=9, frameon=False, loc="upper right")
    ax.grid(axis="y", alpha=0.3, ls=":")
    ax.spines[["top", "right"]].set_visible(False)

style_panel(axes[0], modern_h_years, "#1565C0", "Global Modern Divergences", "A. Global Modern Human Divergence Kinetics")
style_panel(axes[1], squirrel_years, "#C62828", "Red Squirrel Divergences", "B. Red Squirrel Divergence Kinetics", show_fur_trade=True)
style_panel(axes[2], medieval_h_years, "#2E7D32", "Medieval Human Divergences", "C. Medieval Human Divergence Kinetics", post_1100_color="#E65100")

plt.tight_layout(pad=4.0)
plt.savefig(OUT_FIG, dpi=300, bbox_inches="tight")
plt.show()
