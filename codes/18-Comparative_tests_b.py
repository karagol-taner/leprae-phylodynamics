import pandas as pd
import numpy as np
import re
from scipy import stats

# ── CONFIGURATION ─────────────────────────────────────────────
TREE_FILE = "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/lep_MCC_dated.tree"
META_FILE = "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/leprosy_metadata.csv"
REFERENCE_YEAR = 2020

# ── 1. STABLE DATA EXTRACTION (REUSING YOUR PARSER) ───────────
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
        elif token.startswith('[&'):
            h_match = re.search(r'height=([0-9\.E-]+)', token)
            if h_m := h_match:
                if current_node: current_node['height'] = float(h_m.group(1))
        elif token not in '(),;':
            if not token.startswith(':'):
                nid_counter += 1
                tip_node = {'id': nid_counter, 'children': [], 'height': None, 'is_tip': True, 'label': token}
                if stack: stack[-1]['children'].append(tip_node)
                nodes[nid_counter] = tip_node
                current_node = tip_node
    return nodes

def get_descendant_tips(node):
    if node['is_tip']: return {node['label']}
    tips = set()
    for child in node['children']: tips.update(get_descendant_tips(child))
    return tips

# Load Files
with open(TREE_FILE, 'r') as f: content = f.read()
tree_match = re.search(r'tree\s+[^=]+=\s+(.*?;)', content, re.DOTALL | re.IGNORECASE)
nodes_map = parse_nexus_tree(tree_match.group(1))

meta = pd.read_csv(META_FILE)
meta["CleanID"] = meta[meta.columns[0]].astype(str).str.strip()
meta_dict = meta.set_index('CleanID').to_dict('index')

translate_dict = {}
translate_block = re.search(r'Translate\s+(.*?);', content, re.DOTALL | re.IGNORECASE)
if translate_block:
    for line in translate_block.group(1).strip().split('\n'):
        parts = line.strip().split()
        if len(parts) >= 2:
            idx, label = parts[0], parts[1].strip(',').strip("'")
            srr_match = re.search(r'(SRR\d+)', label)
            translate_dict[idx] = srr_match.group(1) if srr_match else label

tip_info = {idx: meta_dict.get(srr, {'Host': 'Unknown'}) for idx, srr in translate_dict.items()}

# Extract Years for All Three Clades
sq_yrs, med_yrs, gl_yrs = [], [], []
for nid, node in nodes_map.items():
    if node['is_tip'] or node['height'] is None: continue
    desc_tips = get_descendant_tips(node)
    total = len(desc_tips)

    sq_idx = {idx for idx, info in tip_info.items() if info['Host'] == 'Red Squirrel'}
    med_idx = {idx for idx, info in tip_info.items() if 'Medieval' in str(info.get('Era','')) and info['Host'] == 'Human'}
    gl_idx = {idx for idx in tip_info if idx not in sq_idx and idx not in med_idx}

    sq_c, med_c, gl_c = len(desc_tips.intersection(sq_idx)), len(desc_tips.intersection(med_idx)), len(desc_tips.intersection(gl_idx))
    year = REFERENCE_YEAR - node['height']
    if sq_c / total > 0.8: sq_yrs.append(year)
    elif med_c / total > 0.8: med_yrs.append(year)
    elif gl_c / total > 0.5: gl_yrs.append(year)

# ── 2. COMPREHENSIVE STATISTICAL SUITE ───────────────────────
def run_final_manuscript_proofs(sq, med, gl):
    # Prepare waiting times (intervals)
    def get_wait(yrs): return np.diff(sorted(yrs))[np.diff(sorted(yrs)) > 1.0]

    w_sq = get_wait(sq)
    w_med = get_wait(med)
    w_gl = get_wait(gl)

    print(f"{'='*70}\n{'EVOLUTIONARY REGIME PROOF: MULTI-HOST COMPARISON':^70}\n{'='*70}")

    # 1. Distribution Tests (K-S Test)
    ks_sq_med = stats.ks_2samp(w_sq, w_med)
    ks_sq_gl = stats.ks_2samp(w_sq, w_gl)

    # 2. Median Tempo Tests (Mann-Whitney)
    u_sq_med = stats.mannwhitneyu(w_sq, w_med)
    u_sq_gl = stats.mannwhitneyu(w_sq, w_gl)

    # 3. Stability/Variance Tests (Levene's)
    lev_sq_med = stats.levene(w_sq, w_med)

    # 4. Effect Size (Cohen's d)
    def cohen_d(x, y): return (np.mean(x) - np.mean(y)) / np.sqrt((np.var(x) + np.var(y)) / 2)
    d_sq_med = cohen_d(w_sq, w_med)
    d_sq_gl = cohen_d(w_sq, w_gl)

    print(f"{'Test Description':<40} | {'p-value':<10} | {'Statistic':<10}")
    print(f"{'-'*70}")
    print(f"{'Identity: Squirrel vs Medieval (K-S)':<40} | {ks_sq_med.pvalue:<10.4f} | {ks_sq_med.statistic:<10.4f}")
    print(f"{'Identity: Squirrel vs Global (K-S)':<40} | {ks_sq_gl.pvalue:<10.4f} | {ks_sq_gl.statistic:<10.4f}")
    print(f"{'Tempo: Squirrel vs Medieval (U-test)':<40} | {u_sq_med.pvalue:<10.4f} | {u_sq_med.statistic:<10.4f}")
    print(f"{'Tempo: Squirrel vs Global (U-test)':<40} | {u_sq_gl.pvalue:<10.4f} | {u_sq_gl.statistic:<10.4f}")
    print(f"{'Stability: Squirrel vs Medieval (Levene)':<40} | {lev_sq_med.pvalue:<10.4f} | {lev_sq_med.statistic:<10.4f}")

    print(f"\n{'EFFECT SIZES (COHEN S D)':-^70}")
    print(f"• Squirrel vs Medieval Human: {d_sq_med:.4f} (Magnitude: {'Large' if abs(d_sq_med)>0.8 else 'Medium'})")
    print(f"• Squirrel vs Global Human:   {d_sq_gl:.4f} (Magnitude: {'Large' if abs(d_sq_gl)>0.8 else 'Medium'})")

    print(f"\n{'SUMMARY INTERPRETATION':-^70}")
    if ks_sq_gl.pvalue < 0.05:
        print("Confirmed: Squirrel reservoir has a distinct regime from Global leprosy.")
    if ks_sq_med.pvalue > 0.05:
        print("Insight: Squirrel reservoir maintains the 'Ancient Regime' of medieval UK leprosy.")
    print(f"{'='*70}")

run_final_manuscript_proofs(sq_yrs, med_yrs, gl_yrs)
