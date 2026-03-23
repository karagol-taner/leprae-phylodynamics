import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import stats
import gc

# --- CONFIGURATION ---
TREE_FILE = "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/lep_MCC_dated.tree"
REFERENCE_YEAR = 2020
START_YEAR = 400

def calculate_gini(x):
    x = np.sort(x)
    n = len(x)
    if n == 0 or np.sum(x) == 0: return np.nan
    index = np.arange(1, n + 1)
    return (np.sum((2 * index - n - 1) * x)) / (n * np.sum(x))

def run_final_manuscript_stats(tree_path):
    # 1. Parsing & Grounding
    tree = Phylo.read(tree_path, "nexus")
    max_dist = max(tree.distance(t) for t in tree.get_terminals())

    # 2. Tip Partitioning
    all_tips = tree.get_terminals()
    sq_tips = [t for t in all_tips if any(s in t.name.lower() for s in ["srr367", "squirrel"])]
    win_tips = [t for t in all_tips if any(w in t.name.lower() for w in ["srr847", "winchester"])]
    global_tips = [t for t in all_tips if t not in sq_tips and t not in win_tips]

    def extract_clade_dates(target_tips, threshold=0.8):
        dates = []
        for node in tree.find_clades(terminal=False):
            node_tips = node.get_terminals()
            if len(set(node_tips).intersection(set(target_tips))) / len(node_tips) >= threshold:
                dist = tree.distance(node)
                year = REFERENCE_YEAR - (max_dist - dist)
                if year >= START_YEAR:
                    dates.append(year)
        return sorted(list(set(dates)))

    # 3. Data Extraction
    sq_dates = extract_clade_dates(sq_tips, 0.8)
    win_dates = extract_clade_dates(win_tips, 0.8)
    gl_dates = extract_clade_dates(global_tips, 0.5)

    def get_metrics(dates):
        if len(dates) < 3: return None
        t = np.sort(dates)
        intervals = np.diff(t)

        # Gamma
        n = len(t) + 1
        T = t - t[0]
        gamma = ((np.sum(T[:-1]) / (n - 2)) - (T[-1] / 2)) / (T[-1] * np.sqrt(1 / (12 * (n - 2))))

        # Saturation (Early vs Late rate ratio)
        mid = len(intervals) // 2
        sat_index = np.mean(intervals[:mid]) / np.mean(intervals[mid:]) if mid > 0 else np.nan

        return {
            "n": len(t), "gamma": gamma, "cv": np.std(intervals) / np.mean(intervals),
            "gini": calculate_gini(intervals), "mean_wait": np.mean(intervals),
            "sat_index": sat_index, "intervals": intervals
        }

    res = {"SQUIRREL": get_metrics(sq_dates), "WINCHESTER": get_metrics(win_dates), "GLOBAL_HUMAN": get_metrics(gl_dates)}

    # --- RESULTS DISPLAY ---
    print(f"\n{'='*75}\n{'FINAL MANUSCRIPT STATISTICAL TABLE':^75}\n{'='*75}")
    print(f"{'Metric':<25} | {'Squirrel':<12} | {'Winchester':<12} | {'Global Human':<12}")
    print("-" * 75)

    metrics = [
        ("Internal Nodes (n)", "n", "{:.0f}"),
        ("Gamma Statistic (γ)", "gamma", "{:.4f}"),
        ("Burstiness (CV)", "cv", "{:.4f}"),
        ("Inequality (Gini)", "gini", "{:.4f}"),
        ("Mean Wait (Years)", "mean_wait", "{:.2f}"),
        ("Saturation Index", "sat_index", "{:.4f}")
    ]

    for label, key, fmt in metrics:
        vals = [fmt.format(res[g][key]) if res[g] else "N/A" for g in ["SQUIRREL", "WINCHESTER", "GLOBAL_HUMAN"]]
        print(f"{label:<25} | {vals[0]:<12} | {vals[1]:<12} | {vals[2]:<12}")

    print(f"\n{'='*75}\n{'COMPARATIVE SIGNIFICANCE TESTS':^75}\n{'='*75}")
    if res["SQUIRREL"] and res["GLOBAL_HUMAN"]:
        u, p = stats.mannwhitneyu(res["SQUIRREL"]["intervals"], res["GLOBAL_HUMAN"]["intervals"])
        print(f"• Tempo Difference (Sq vs Global): p = {p:.4f}")
        print(f"• Diversification Ratio (Sq/Gl):   {res['GLOBAL_HUMAN']['mean_wait']/res['SQUIRREL']['mean_wait']:.2f}x")

    if res["SQUIRREL"]:
        status = "EXPANDING" if res["SQUIRREL"]['sat_index'] > 1.2 else "STABLE" if 0.8 <= res["SQUIRREL"]['sat_index'] <= 1.2 else "CONTRACTING"
        print(f"• Squirrel Reservoir Status:      {status}")
    print(f"{'='*75}")

run_final_manuscript_stats(TREE_FILE)
