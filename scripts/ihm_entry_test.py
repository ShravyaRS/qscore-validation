"""
Q-score computation on a PDB-IHM integrative structure.
Demonstrates that the pure Python Q-score works on IHMValidation's
target data type: integrative models with associated cryo-EM maps.

Entry: PDBDEV_00000141 / EMD-14774 (PTX3 Pentraxin, 2.5 Å)
"""

import numpy as np
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

def main():
    cif_path = 'test_data/PDBDEV_00000141/PDBDEV_00000141.cif'
    map_path = 'test_data/PDBDEV_00000141/emd_14774.map.gz'
    
    print(f"\n{'='*70}")
    print(f"Q-SCORE ON PDB-IHM ENTRY: PDBDEV_00000141 / EMD-14774")
    print(f"PTX3 Pentraxin Domain (2.5 Å, integrative model)")
    print(f"{'='*70}")
    
    # Load structure info
    from qscore.pdb_utils import get_protein_from_file_path
    from qscore.mrc_utils import load_mrc
    from qscore.q_score import calculate_q_score
    
    print("\n[1/4] Loading structure...")
    protein = get_protein_from_file_path(cif_path)
    n_res = protein.aatype.shape[0]
    n_atoms = int(protein.atom_mask.sum())
    n_chains = len(np.unique(protein.chain_index))
    print(f"  → {n_res} residues, {n_atoms} atoms, {n_chains} chains")
    
    print("\n[2/4] Loading map...")
    mrc = load_mrc(map_path)
    print(f"  → Grid: {mrc.grid.shape}")
    print(f"  → Voxel: {mrc.voxel_size}")
    print(f"  → Map stats: mean={mrc.grid.mean():.4f}, std={mrc.grid.std():.4f}, max={mrc.grid.max():.4f}")
    
    print("\n[3/4] Computing Q-scores (σ=0.4)...")
    atoms_flat = protein.atom_positions.reshape(-1, 3)
    mask_flat = protein.atom_mask.reshape(-1)
    valid_atoms = atoms_flat[mask_flat > 0]
    
    atom_qs = calculate_q_score(valid_atoms, mrc, ref_gaussian_width=0.4)
    
    # Per-residue aggregation
    residue_qs = {}
    chain_qs = {}
    atom_idx = 0
    
    for i in range(n_res):
        n_atoms_res = int(protein.atom_mask[i].sum())
        chain_idx = int(protein.chain_index[i])
        resnum = int(protein.residue_index[i])
        
        if n_atoms_res > 0 and atom_idx + n_atoms_res <= len(atom_qs):
            res_q = atom_qs[atom_idx:atom_idx + n_atoms_res]
            valid = res_q[~np.isnan(res_q)]
            if len(valid) > 0:
                mean_q = float(np.mean(valid))
                
                # Get chain ID (chain_id is per-chain, chain_index maps residue→chain)
                chain_id = str(protein.chain_id[chain_idx])
                
                residue_qs[(chain_id, resnum)] = mean_q
                
                if chain_id not in chain_qs:
                    chain_qs[chain_id] = []
                chain_qs[chain_id].append(mean_q)
            
            atom_idx += n_atoms_res
    
    # Overall stats
    all_res_qs = list(residue_qs.values())
    overall_avg = np.mean(all_res_qs)
    
    print(f"\n{'='*70}")
    print(f"RESULTS")
    print(f"{'='*70}")
    print(f"\n  Overall average Q-score: {overall_avg:.4f}")
    print(f"  Residues computed:       {len(residue_qs)}")
    print(f"  Q-score range:           [{min(all_res_qs):.3f}, {max(all_res_qs):.3f}]")
    print(f"  Median:                  {np.median(all_res_qs):.4f}")
    print(f"  Std:                     {np.std(all_res_qs):.4f}")
    
    print(f"\n  Per-chain Q-scores:")
    print(f"  {'Chain':<8} {'N_res':<8} {'Mean Q':<10} {'Std':<10} {'Min':<10} {'Max':<10}")
    print(f"  {'-'*56}")
    for chain_id in sorted(chain_qs.keys()):
        qs = chain_qs[chain_id]
        print(f"  {chain_id:<8} {len(qs):<8} {np.mean(qs):<10.4f} {np.std(qs):<10.4f} "
              f"{min(qs):<10.4f} {max(qs):<10.4f}")
    
    # Save results
    results = {
        'entry': 'PDBDEV_00000141',
        'emdb': 'EMD-14774',
        'pdb': '9A24',
        'resolution': 2.5,
        'description': 'PTX3 Pentraxin Domain (integrative: cryo-EM + AlphaFold)',
        'n_residues': len(residue_qs),
        'n_atoms': n_atoms,
        'n_chains': n_chains,
        'overall_qscore': round(overall_avg, 4),
        'qscore_median': round(float(np.median(all_res_qs)), 4),
        'qscore_std': round(float(np.std(all_res_qs)), 4),
        'qscore_range': [round(min(all_res_qs), 4), round(max(all_res_qs), 4)],
        'per_chain': {
            chain_id: {
                'n_residues': len(qs),
                'mean_qscore': round(float(np.mean(qs)), 4),
                'std_qscore': round(float(np.std(qs)), 4),
            }
            for chain_id, qs in sorted(chain_qs.items())
        }
    }
    
    json_path = 'results/ihm_pdbdev_00000141.json'
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved: {json_path}")
    
    # Generate figure: Q-score distribution + per-chain box plot
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    
    fig = plt.figure(figsize=(14, 5))
    gs = GridSpec(1, 3, width_ratios=[1, 1.2, 1], wspace=0.35)
    
    # Panel A: Q-score distribution
    ax1 = fig.add_subplot(gs[0])
    ax1.hist(all_res_qs, bins=40, color='#2563EB', alpha=0.7, edgecolor='white', linewidth=0.3)
    ax1.axvline(overall_avg, color='red', lw=1.5, ls='--', label=f'Mean = {overall_avg:.3f}')
    ax1.set_xlabel('Q-score', fontsize=10)
    ax1.set_ylabel('Count', fontsize=10)
    ax1.set_title('A. Q-score distribution', fontsize=11, fontweight='bold')
    ax1.legend(fontsize=8)
    
    # Panel B: Per-chain box plot
    ax2 = fig.add_subplot(gs[1])
    chain_ids_sorted = sorted(chain_qs.keys())
    chain_data = [chain_qs[c] for c in chain_ids_sorted]
    bp = ax2.boxplot(chain_data, tick_labels=chain_ids_sorted, patch_artist=True,
                     medianprops=dict(color='red', lw=1.5))
    for patch in bp['boxes']:
        patch.set_facecolor('#2563EB')
        patch.set_alpha(0.6)
    ax2.set_xlabel('Chain', fontsize=10)
    ax2.set_ylabel('Q-score', fontsize=10)
    ax2.set_title('B. Per-chain Q-scores', fontsize=11, fontweight='bold')
    ax2.tick_params(axis='x', rotation=45)
    
    # Panel C: Q-score along sequence (first chain)
    ax3 = fig.add_subplot(gs[2])
    first_chain = chain_ids_sorted[0]
    chain_residues = [(k[1], v) for k, v in residue_qs.items() if k[0] == first_chain]
    chain_residues.sort()
    if chain_residues:
        resnums, qs_vals = zip(*chain_residues)
        ax3.plot(resnums, qs_vals, 'o', markersize=3, alpha=0.5, color='#2563EB')
        # Rolling average
        window = 10
        if len(qs_vals) > window:
            rolling = np.convolve(qs_vals, np.ones(window)/window, mode='valid')
            rolling_x = resnums[window//2:window//2+len(rolling)]
            ax3.plot(rolling_x, rolling, 'r-', lw=1.5, label=f'{window}-res moving avg')
        ax3.set_xlabel(f'Residue number (Chain {first_chain})', fontsize=10)
        ax3.set_ylabel('Q-score', fontsize=10)
        ax3.set_title(f'C. Q-score along sequence', fontsize=11, fontweight='bold')
        ax3.legend(fontsize=8)
    
    plt.suptitle('PDB-IHM Entry: PDBDEV_00000141 / EMD-14774 (PTX3, 2.5 Å)',
                 fontsize=13, fontweight='bold', y=1.02)
    
    fig_path = 'results/ihm_pdbdev_00000141.png'
    plt.savefig(fig_path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fig_path}")
    
    print(f"\n{'='*70}")
    print("DONE — Pure Python Q-score works on PDB-IHM integrative structures")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
