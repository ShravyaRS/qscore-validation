"""
Per-residue Q-score comparison: wwPDB reference (MapQ) vs jamaliki/qscore
Validates agreement at the individual residue level, not just global averages.
"""

import gzip
import xml.etree.ElementTree as ET
import numpy as np
import mrcfile
import os
import sys
import json
from pathlib import Path

# ---------------------------------------------------------------------------
# 1. Parse per-residue reference Q-scores from wwPDB validation XML
# ---------------------------------------------------------------------------
def parse_wwpdb_qscores(xml_gz_path):
    """Extract per-residue Q-scores from wwPDB validation XML."""
    with gzip.open(xml_gz_path, 'rt') as f:
        tree = ET.parse(f)
    root = tree.getroot()
    
    residues = {}
    for elem in root.iter('ModelledSubgroup'):
        attrib = elem.attrib
        if 'Q_score' in attrib and 'chain' in attrib and 'resnum' in attrib:
            chain = attrib['chain']
            resnum = int(attrib['resnum'])
            resname = attrib.get('resname', 'UNK')
            qscore = float(attrib['Q_score'])
            model = attrib.get('model', '1')
            if model == '1':  # Only model 1
                key = (chain, resnum)
                residues[key] = {
                    'chain': chain,
                    'resnum': resnum,
                    'resname': resname,
                    'qscore_mapq': qscore
                }
    return residues

# ---------------------------------------------------------------------------
# 2. Compute per-residue Q-scores using jamaliki/qscore
# ---------------------------------------------------------------------------
def compute_per_residue_qscores(cif_path, map_path):
    """Compute per-residue Q-scores using qscore package."""
    from qscore.q_score import calculate_q_score
    from qscore.pdb_utils import get_protein_from_file_path
    from qscore.mrc_utils import load_mrc
    
    # Load structure and map using qscore's own loaders
    protein = get_protein_from_file_path(cif_path)
    mrc = load_mrc(map_path)
    
    # Flatten atom positions and mask (only real atoms)
    atoms_flat = protein.atom_positions.reshape(-1, 3)
    mask_flat = protein.atom_mask.reshape(-1)
    valid_atoms = atoms_flat[mask_flat > 0]
    
    # Compute per-atom Q-scores (sigma=0.4 → ref_gaussian_width=0.4)
    atom_qs = calculate_q_score(valid_atoms, mrc, ref_gaussian_width=0.4)
    
    # Map back to residues using protein structure
    residue_scores = {}
    atom_idx = 0
    n_res = protein.aatype.shape[0]
    
    for i in range(n_res):
        chain_idx = int(protein.chain_index[i])
        resnum = int(protein.residue_index[i])
        n_atoms_res = int(protein.atom_mask[i].sum())
        
        if n_atoms_res > 0 and atom_idx + n_atoms_res <= len(atom_qs):
            res_q = atom_qs[atom_idx:atom_idx + n_atoms_res]
            valid = res_q[~np.isnan(res_q)]
            if len(valid) > 0:
                # Get chain ID
                chain_id = chr(65 + chain_idx) if chain_idx < 26 else str(chain_idx)
                if hasattr(protein, 'chain_id') and len(protein.chain_id) > i:
                    cid = protein.chain_id[i]
                    if isinstance(cid, (bytes, np.bytes_)):
                        chain_id = cid.decode('utf-8') if isinstance(cid, bytes) else str(cid)
                    elif isinstance(cid, str):
                        chain_id = cid
                    else:
                        chain_id = str(cid)
                
                res_key = (chain_id, resnum)
                residue_scores[res_key] = {
                    'chain': chain_id,
                    'resnum': resnum,
                    'resname': '',
                    'qscore_python': float(np.mean(valid)),
                    'n_atoms': len(valid)
                }
            atom_idx += n_atoms_res
    
    return residue_scores

# ---------------------------------------------------------------------------
# 3. Main comparison
# ---------------------------------------------------------------------------
def main():
    # Paths
    entry = '5a1a'
    emdb = '2984'
    xml_gz = '/tmp/5a1a_val.xml.gz'
    
    test_dir = Path('test_data') / f'{entry}'
    cif_path = str(test_dir / f'{entry}.cif')
    map_path = str(test_dir / f'emd_{emdb}.map.gz')
    
    # Check files exist
    if not os.path.exists(xml_gz):
        print(f"Downloading validation XML...")
        os.system(f'wget -q "https://ftp.ebi.ac.uk/pub/databases/pdb/validation_reports/a1/{entry}/{entry}_validation.xml.gz" -O {xml_gz}')
    
    if not os.path.exists(cif_path):
        print(f"Downloading structure files...")
        os.makedirs(test_dir, exist_ok=True)
        os.system(f'wget -q "https://files.rcsb.org/download/{entry}.cif" -O {cif_path}')
    
    if not os.path.exists(map_path):
        os.system(f'wget -q "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{emdb}/map/emd_{emdb}.map.gz" -O {map_path}')
    
    print(f"\n{'='*70}")
    print(f"PER-RESIDUE Q-SCORE COMPARISON: {entry.upper()} / EMD-{emdb}")
    print(f"{'='*70}")
    
    # Parse reference
    print("\n[1/3] Parsing wwPDB reference Q-scores...")
    ref_scores = parse_wwpdb_qscores(xml_gz)
    print(f"  → {len(ref_scores)} residues with reference Q-scores")
    
    # Compute qscore
    print("\n[2/3] Computing Q-scores with jamaliki/qscore (σ=0.4)...")
    py_scores = compute_per_residue_qscores(cif_path, map_path)
    print(f"  → {len(py_scores)} residues computed")
    
    # Match residues
    print("\n[3/3] Matching and comparing...")
    matched_keys = set(ref_scores.keys()) & set(py_scores.keys())
    print(f"  → {len(matched_keys)} residues matched")
    
    mapq_vals = []
    python_vals = []
    rows = []
    
    for key in sorted(matched_keys):
        ref = ref_scores[key]
        py = py_scores[key]
        mapq_vals.append(ref['qscore_mapq'])
        python_vals.append(py['qscore_python'])
        rows.append({
            'chain': key[0],
            'resnum': key[1],
            'resname': py['resname'],
            'qscore_mapq': ref['qscore_mapq'],
            'qscore_python': py['qscore_python'],
            'diff': py['qscore_python'] - ref['qscore_mapq'],
            'n_atoms': py['n_atoms']
        })
    
    mapq_arr = np.array(mapq_vals)
    py_arr = np.array(python_vals)
    diff_arr = py_arr - mapq_arr
    
    # Statistics
    from scipy import stats
    pearson_r, pearson_p = stats.pearsonr(mapq_arr, py_arr)
    spearman_r, spearman_p = stats.spearmanr(mapq_arr, py_arr)
    slope, intercept, _, _, _ = stats.linregress(mapq_arr, py_arr)
    rmsd = np.sqrt(np.mean(diff_arr**2))
    
    print(f"\n{'='*70}")
    print(f"RESULTS: Per-residue comparison ({len(matched_keys)} residues)")
    print(f"{'='*70}")
    print(f"  Pearson r     = {pearson_r:.4f}  (p = {pearson_p:.2e})")
    print(f"  Spearman ρ    = {spearman_r:.4f}  (p = {spearman_p:.2e})")
    print(f"  Mean offset   = {np.mean(diff_arr):+.4f} ± {np.std(diff_arr):.4f}")
    print(f"  Median offset = {np.median(diff_arr):+.4f}")
    print(f"  RMSD          = {rmsd:.4f}")
    print(f"  Max |diff|    = {np.max(np.abs(diff_arr)):.4f}")
    print(f"  Linear fit    = {slope:.4f} × MapQ + {intercept:.4f}")
    print(f"  R²            = {pearson_r**2:.4f}")
    
    # Save CSV
    csv_path = 'results/per_residue_5a1a.csv'
    with open(csv_path, 'w') as f:
        f.write('chain,resnum,resname,qscore_mapq,qscore_python,diff,n_atoms\n')
        for row in rows:
            f.write(f"{row['chain']},{row['resnum']},{row['resname']},"
                    f"{row['qscore_mapq']:.4f},{row['qscore_python']:.4f},"
                    f"{row['diff']:.4f},{row['n_atoms']}\n")
    print(f"\n  Saved: {csv_path}")
    
    # ---------------------------------------------------------------------------
    # 4. Generate publication-quality figure
    # ---------------------------------------------------------------------------
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    
    fig = plt.figure(figsize=(16, 5))
    gs = GridSpec(1, 3, width_ratios=[1, 1, 1], wspace=0.35)
    
    # Panel A: Scatter plot
    ax1 = fig.add_subplot(gs[0])
    ax1.scatter(mapq_arr, py_arr, s=4, alpha=0.3, c='#2563EB', edgecolors='none')
    lims = [min(mapq_arr.min(), py_arr.min()) - 0.05, 
            max(mapq_arr.max(), py_arr.max()) + 0.05]
    ax1.plot(lims, lims, 'k--', lw=0.8, alpha=0.5, label='y = x')
    fit_x = np.linspace(lims[0], lims[1], 100)
    ax1.plot(fit_x, slope * fit_x + intercept, 'r-', lw=1.2, 
             label=f'Fit: y = {slope:.3f}x + {intercept:.3f}')
    ax1.set_xlabel('MapQ Q-score (wwPDB reference)', fontsize=10)
    ax1.set_ylabel('Python Q-score (jamaliki/qscore)', fontsize=10)
    ax1.set_title(f'A. Per-residue correlation (n={len(matched_keys)})', fontsize=11, fontweight='bold')
    ax1.legend(fontsize=8, loc='upper left')
    ax1.set_xlim(lims)
    ax1.set_ylim(lims)
    ax1.set_aspect('equal')
    ax1.text(0.95, 0.05, f'r = {pearson_r:.4f}\nρ = {spearman_r:.4f}', 
             transform=ax1.transAxes, ha='right', va='bottom', fontsize=9,
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    # Panel B: Histogram of differences
    ax2 = fig.add_subplot(gs[1])
    ax2.hist(diff_arr, bins=60, color='#2563EB', alpha=0.7, edgecolor='white', linewidth=0.3)
    ax2.axvline(np.mean(diff_arr), color='red', lw=1.5, ls='--', 
                label=f'Mean = {np.mean(diff_arr):+.4f}')
    ax2.axvline(0, color='black', lw=0.8, ls='-', alpha=0.5)
    ax2.set_xlabel('Q-score difference (Python − MapQ)', fontsize=10)
    ax2.set_ylabel('Count', fontsize=10)
    ax2.set_title('B. Distribution of residue-level differences', fontsize=11, fontweight='bold')
    ax2.legend(fontsize=8)
    
    # Panel C: Difference by residue position (chain A only)
    ax3 = fig.add_subplot(gs[2])
    chain_a = [(r['resnum'], r['diff']) for r in rows if r['chain'] == 'A']
    if chain_a:
        resnums, diffs = zip(*sorted(chain_a))
        ax3.scatter(resnums, diffs, s=6, alpha=0.5, c='#2563EB', edgecolors='none')
        # Rolling average
        window = 20
        if len(diffs) > window:
            rolling = np.convolve(diffs, np.ones(window)/window, mode='valid')
            rolling_x = resnums[window//2:window//2+len(rolling)]
            ax3.plot(rolling_x, rolling, 'r-', lw=1.5, label=f'{window}-residue moving avg')
        ax3.axhline(0, color='black', lw=0.8, ls='-', alpha=0.5)
        ax3.set_xlabel('Residue number (Chain A)', fontsize=10)
        ax3.set_ylabel('Q-score difference', fontsize=10)
        ax3.set_title('C. Offset along sequence (Chain A)', fontsize=11, fontweight='bold')
        ax3.legend(fontsize=8)
    
    plt.suptitle(f'Per-residue Q-score validation: {entry.upper()} / EMD-{emdb} (2.2 Å)',
                 fontsize=13, fontweight='bold', y=1.02)
    
    fig_path = 'results/per_residue_correlation.png'
    plt.savefig(fig_path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fig_path}")
    
    # Summary JSON
    summary = {
        'entry': entry.upper(),
        'emdb': f'EMD-{emdb}',
        'resolution': 2.2,
        'n_residues_matched': len(matched_keys),
        'n_residues_ref': len(ref_scores),
        'n_residues_computed': len(py_scores),
        'pearson_r': round(pearson_r, 4),
        'spearman_rho': round(spearman_r, 4),
        'mean_offset': round(float(np.mean(diff_arr)), 4),
        'std_offset': round(float(np.std(diff_arr)), 4),
        'median_offset': round(float(np.median(diff_arr)), 4),
        'rmsd': round(rmsd, 4),
        'max_abs_diff': round(float(np.max(np.abs(diff_arr))), 4),
        'linear_fit': {'slope': round(slope, 4), 'intercept': round(intercept, 4)},
        'r_squared': round(pearson_r**2, 4)
    }
    
    json_path = 'results/per_residue_summary.json'
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"  Saved: {json_path}")
    
    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
