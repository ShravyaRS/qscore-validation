#!/usr/bin/env python3
"""
Q-score Validation: Compare jamaliki/qscore (pure Python) vs MapQ/Chimera (EMDB reference)

Downloads PDB structures and EMDB maps, computes Q-scores using the pure Python
reimplementation, and compares against EMDB reference values.

References:
    [1] Pintilie et al. (2020). Nat Methods 17, 328-334.
    [2] Pintilie et al. (2025). Acta Cryst. D81.
"""

import os
import csv
import subprocess
import urllib.request
import json
import numpy as np
from scipy.stats import spearmanr, pearsonr


def get_test_entries(n=30):
    """Fetch EMDB entries with Q-scores and their PDB codes."""
    url = ("https://www.ebi.ac.uk/emdb/api/search/"
           "resolution:[1.5 TO 5.0] AND average_qscore_value:[0.1 TO *]"
           f"?rows={n*8}&wt=csv&fl=emdb_id,resolution,average_qscore_value")
    response = urllib.request.urlopen(url, timeout=30)
    lines = response.read().decode().strip().split("\n")
    
    entries = []
    reader = csv.DictReader(lines)
    for row in reader:
        q = row["average_qscore_value"]
        if "," in q:
            continue
        entries.append({
            "emdb": row["emdb_id"],
            "res": float(row["resolution"]),
            "ref_q": float(q),
        })
    
    entries.sort(key=lambda x: x["res"])
    step = max(1, len(entries) // n)
    selected = entries[::step][:n]
    
    # Get PDB codes
    for e in selected:
        try:
            api_url = f"https://www.ebi.ac.uk/emdb/api/entry/{e['emdb']}"
            req = urllib.request.urlopen(api_url, timeout=10)
            data = json.loads(req.read())
            pdb_list = data.get("crossreferences", {}).get("pdb_list", {}).get("pdb_reference", [])
            if isinstance(pdb_list, dict):
                pdb_list = [pdb_list]
            e["pdb"] = pdb_list[0].get("pdb_id", "") if pdb_list else ""
        except Exception:
            e["pdb"] = ""
    
    return [e for e in selected if e["pdb"]]


def download_files(emdb_id, pdb_id, data_dir="test_data"):
    """Download PDB structure and EMDB map."""
    os.makedirs(data_dir, exist_ok=True)
    num = emdb_id.replace("EMD-", "")
    cif = os.path.join(data_dir, f"{pdb_id.upper()}.cif")
    mapgz = os.path.join(data_dir, f"emd_{num}.map.gz")
    mapf = os.path.join(data_dir, f"emd_{num}.map")

    if not os.path.exists(cif) or os.path.getsize(cif) == 0:
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
        subprocess.run(["wget", "-q", url, "-O", cif], check=True)

    if not os.path.exists(mapf):
        if not os.path.exists(mapgz) or os.path.getsize(mapgz) == 0:
            url = f"https://ftp.ebi.ac.uk/pub/databases/emdb/structures/{emdb_id}/map/emd_{num}.map.gz"
            subprocess.run(["wget", "-q", url, "-O", mapgz], check=True)
        subprocess.run(["gunzip", "-f", mapgz], check=True)

    return cif, mapf


def check_alignment(cif_path, map_path):
    """Verify atoms fall within map bounds."""
    import mrcfile
    from Bio.PDB import MMCIFParser

    mrc = mrcfile.open(map_path, "r")
    voxel = float(mrc.voxel_size.x)
    origin = np.array([mrc.header.origin.x, mrc.header.origin.y, mrc.header.origin.z])
    nxs = np.array([int(mrc.header.nxstart), int(mrc.header.nystart), int(mrc.header.nzstart)])
    real_origin = origin + nxs * voxel
    grid_max = real_origin + np.array(mrc.data.shape[::-1]) * voxel
    mrc.close()

    parser = MMCIFParser(QUIET=True)
    s = parser.get_structure("t", cif_path)
    atoms = np.array([a.get_vector().get_array() for a in s.get_atoms()])
    return np.all(atoms >= real_origin - 2) and np.all(atoms <= grid_max + 2)


def compute_qscore(cif_path, map_path, sigma=0.4):
    """Compute mean Q-score using jamaliki/qscore."""
    from qscore.q_score import calculate_q_score
    from qscore.mrc_utils import load_mrc
    from qscore.pdb_utils import get_protein_from_file_path

    prot = get_protein_from_file_path(cif_path)
    mrc_map = load_mrc(map_path, False)
    atoms = prot.atom_positions[prot.atom_mask.astype(bool)]
    np.random.seed(42)
    q = calculate_q_score(atoms, mrc_map, ref_gaussian_width=sigma)
    return float(np.mean(q))


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--entries", type=int, default=25, help="Number of entries to test")
    args = parser.parse_args()

    print("Fetching EMDB entries with Q-scores...")
    entries = get_test_entries(args.entries)
    print(f"Found {len(entries)} entries with PDB codes")

    results = []
    for e in entries:
        emdb, pdb, res, ref_q = e["emdb"], e["pdb"].upper(), e["res"], e["ref_q"]
        print(f"\n{emdb}/{pdb} (res={res:.1f} Å):")

        try:
            cif, mapf = download_files(emdb, pdb)
            if not check_alignment(cif, mapf):
                print("  SKIP: atoms outside map bounds")
                continue
            our_q = compute_qscore(cif, mapf, sigma=0.4)
            diff = our_q - ref_q
            print(f"  MapQ={ref_q:.3f}, qscore={our_q:.4f}, diff={diff:+.4f}")
            results.append((emdb, pdb, res, ref_q, our_q, diff))
        except Exception as ex:
            print(f"  ERROR: {ex}")

    # Save CSV
    os.makedirs("results", exist_ok=True)
    with open("results/qscore_comparison_full.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["emdb_id", "pdb_id", "resolution", "mapq_qscore",
                     "python_qscore", "difference", "pct_difference"])
        for r in sorted(results, key=lambda x: x[2]):
            pct = (r[5] / r[3]) * 100
            w.writerow([r[0], r[1], r[2], r[3], f"{r[4]:.4f}", f"{r[5]:.4f}", f"{pct:.1f}"])

    # Stats
    ref = np.array([r[3] for r in results])
    ours = np.array([r[4] for r in results])
    diffs = ours - ref
    pr, _ = pearsonr(ref, ours)
    sr, _ = spearmanr(ref, ours)
    print(f"\n{'='*60}")
    print(f"RESULTS (n={len(results)})")
    print(f"  Pearson r:    {pr:.4f}")
    print(f"  Spearman r:   {sr:.4f}")
    print(f"  Mean offset:  {np.mean(diffs):+.4f} ± {np.std(diffs):.4f}")
    print(f"  Max |diff|:   {np.max(np.abs(diffs)):.4f}")


if __name__ == "__main__":
    main()
