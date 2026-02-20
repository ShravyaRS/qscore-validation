#!/usr/bin/env python3
"""
Analyze Q-score comparison results and generate plots.

References:
    [1] Pintilie et al. (2020). Nat Methods 17, 328-334.
    [2] Pintilie et al. (2025). Acta Cryst. D81.
    [3] Lawson et al. (2021). Nat Methods 18, 156-164.
"""

import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

def load_results(csvfile):
    results = []
    with open(csvfile) as f:
        for row in csv.DictReader(f):
            results.append({
                "emdb": row["emdb_id"],
                "pdb": row["pdb_id"],
                "res": float(row["resolution"]),
                "mapq": float(row["mapq_qscore"]),
                "ours": float(row["python_qscore"]),
                "diff": float(row["difference"]),
            })
    return results

def identify_outliers(results, threshold=0.1):
    outliers = [r for r in results if abs(r["diff"]) > threshold]
    clean = [r for r in results if abs(r["diff"]) <= threshold]
    return clean, outliers

def compute_stats(results, label=""):
    ref = np.array([r["mapq"] for r in results])
    ours = np.array([r["ours"] for r in results])
    diffs = ours - ref
    pr, pp = pearsonr(ref, ours)
    sr, sp = spearmanr(ref, ours)
    slope, intercept = np.polyfit(ref, ours, 1)
    return {
        "n": len(results),
        "mean_offset": np.mean(diffs),
        "std_offset": np.std(diffs),
        "max_abs_diff": np.max(np.abs(diffs)),
        "pearson_r": pr, "pearson_p": pp,
        "spearman_r": sr, "spearman_p": sp,
        "slope": slope, "intercept": intercept,
        "label": label,
    }

def plot_correlation(results, outliers, outfile):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    clean_ref = np.array([r["mapq"] for r in results])
    clean_ours = np.array([r["ours"] for r in results])
    out_ref = np.array([r["mapq"] for r in outliers]) if outliers else np.array([])
    out_ours = np.array([r["ours"] for r in outliers]) if outliers else np.array([])

    # Panel A: Correlation scatter
    ax = axes[0]
    ax.scatter(clean_ref, clean_ours, c="steelblue", s=60, edgecolors="k", linewidth=0.5, label=f"Valid (n={len(results)})", zorder=3)
    if len(outliers) > 0:
        ax.scatter(out_ref, out_ours, c="red", s=60, marker="x", linewidth=2, label=f"Outlier (n={len(outliers)})", zorder=3)
    lims = [min(clean_ref.min(), clean_ours.min()) - 0.05, max(clean_ref.max(), clean_ours.max()) + 0.05]
    ax.plot(lims, lims, "k--", alpha=0.4, label="y=x")
    slope, intercept = np.polyfit(clean_ref, clean_ours, 1)
    x_fit = np.linspace(lims[0], lims[1], 100)
    ax.plot(x_fit, slope * x_fit + intercept, "r-", alpha=0.7, label=f"Fit: y={slope:.3f}x+{intercept:.3f}")
    pr, _ = pearsonr(clean_ref, clean_ours)
    sr, _ = spearmanr(clean_ref, clean_ours)
    ax.set_xlabel("MapQ Q-score (EMDB reference)", fontsize=12)
    ax.set_ylabel("Python qscore", fontsize=12)
    ax.set_title(f"A. Correlation (Pearson r={pr:.3f}, Spearman ρ={sr:.3f})", fontsize=12)
    ax.legend(fontsize=9)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    # Panel B: Difference vs resolution
    ax = axes[1]
    res_clean = np.array([r["res"] for r in results])
    diffs_clean = clean_ours - clean_ref
    ax.scatter(res_clean, diffs_clean, c="steelblue", s=60, edgecolors="k", linewidth=0.5, zorder=3)
    if len(outliers) > 0:
        res_out = np.array([r["res"] for r in outliers])
        diffs_out = out_ours - out_ref
        ax.scatter(res_out, diffs_out, c="red", s=60, marker="x", linewidth=2, zorder=3)
    ax.axhline(y=np.mean(diffs_clean), color="r", linestyle="--", alpha=0.7, label=f"Mean offset: {np.mean(diffs_clean):+.3f}")
    ax.axhline(y=0, color="k", linestyle="-", alpha=0.3)
    ax.fill_between([1, 5], np.mean(diffs_clean) - np.std(diffs_clean),
                    np.mean(diffs_clean) + np.std(diffs_clean), alpha=0.1, color="red")
    ax.set_xlabel("Resolution (Å)", fontsize=12)
    ax.set_ylabel("Q-score difference (qscore − MapQ)", fontsize=12)
    ax.set_title("B. Offset vs Resolution", fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel C: Bland-Altman plot
    ax = axes[2]
    mean_vals = (clean_ref + clean_ours) / 2
    ax.scatter(mean_vals, diffs_clean, c="steelblue", s=60, edgecolors="k", linewidth=0.5, zorder=3)
    mean_d = np.mean(diffs_clean)
    std_d = np.std(diffs_clean)
    ax.axhline(y=mean_d, color="r", linestyle="--", label=f"Mean: {mean_d:+.4f}")
    ax.axhline(y=mean_d + 1.96*std_d, color="gray", linestyle=":", label=f"+1.96σ: {mean_d+1.96*std_d:+.4f}")
    ax.axhline(y=mean_d - 1.96*std_d, color="gray", linestyle=":", label=f"−1.96σ: {mean_d-1.96*std_d:+.4f}")
    ax.set_xlabel("Mean Q-score (MapQ + qscore) / 2", fontsize=12)
    ax.set_ylabel("Difference (qscore − MapQ)", fontsize=12)
    ax.set_title("C. Bland-Altman Plot", fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"  Plot saved: {outfile}")

def main():
    results = load_results("results/qscore_comparison_full.csv")
    clean, outliers = identify_outliers(results)

    print("=" * 70)
    print("Q-SCORE VALIDATION ANALYSIS")
    print("jamaliki/qscore vs MapQ/Chimera (EMDB reference)")
    print("=" * 70)

    # All data
    stats_all = compute_stats(results, "All entries")
    print(f"\nALL ENTRIES (n={stats_all['n']}):")
    print(f"  Pearson r:   {stats_all['pearson_r']:.4f}")
    print(f"  Spearman r:  {stats_all['spearman_r']:.4f}")
    print(f"  Mean offset: {stats_all['mean_offset']:+.4f} ± {stats_all['std_offset']:.4f}")
    print(f"  Max |diff|:  {stats_all['max_abs_diff']:.4f}")

    # Without outliers
    stats_clean = compute_stats(clean, "Without outliers")
    print(f"\nWITHOUT OUTLIERS (n={stats_clean['n']}):")
    print(f"  Pearson r:   {stats_clean['pearson_r']:.4f}")
    print(f"  Spearman r:  {stats_clean['spearman_r']:.4f}")
    print(f"  Mean offset: {stats_clean['mean_offset']:+.4f} ± {stats_clean['std_offset']:.4f}")
    print(f"  Max |diff|:  {stats_clean['max_abs_diff']:.4f}")
    print(f"  Linear fit:  qscore = {stats_clean['slope']:.3f} × MapQ + {stats_clean['intercept']:.3f}")

    if outliers:
        print(f"\nOUTLIERS ({len(outliers)}):")
        for o in outliers:
            print(f"  {o['emdb']}/{o['pdb']}: res={o['res']}Å, MapQ={o['mapq']:.3f}, qscore={o['ours']:.4f}, diff={o['diff']:+.4f}")
        print(f"  Cause: extremely low map contrast (map std ~100x below normal)")

    # Generate plots
    plot_correlation(clean, outliers, "results/qscore_correlation.png")

    # Print table
    print(f"\n{'Entry':<25} {'Res(Å)':<8} {'MapQ':<8} {'qscore':<8} {'Diff':<10} {'%Diff':<8}")
    print("-" * 70)
    for r in sorted(results, key=lambda x: x["res"]):
        flag = " *" if r in outliers else ""
        pct = (r["diff"] / r["mapq"]) * 100
        print(f"{r['emdb']}/{r['pdb']:<10} {r['res']:<8.1f} {r['mapq']:<8.3f} {r['ours']:<8.4f} {r['diff']:+.4f}    {pct:+.1f}%{flag}")

if __name__ == "__main__":
    main()
