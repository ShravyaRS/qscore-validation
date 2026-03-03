# Complete Analysis: jamaliki/qscore vs MapQ

## Result

When all implementation differences are addressed, our pure-Python Q-score
reproduces MapQ to within ±0.1% (±0.0004 absolute).

| Entry | Our Q (protein) | MapQ | Difference |
|-------|----------------|------|------------|
| EMD-2984/5A1A (2.2Å) | 0.6154 | 0.615 | +0.0004 (+0.1%) |
| EMD-72359/9XZK (3.9Å) | 0.3764 | 0.376 | +0.0004 (+0.1%) |

## Six differences identified

### 1. Epsilon stabilization (critical)
- jamaliki: eps=1e-6 in Pearson denominator
- MapQ: no epsilon
- Fix: use denom_sq > 0 guard instead of additive epsilon

### 2. Sphere sampling
- jamaliki: random (np.random.randn), non-deterministic
- MapQ: Saff-Kuijlaars spiral, deterministic

### 3. Neighbor exclusion
- jamaliki: point closest to originating atom (KDTree nearest)
- MapQ: nearest neighbor distance > 0.9 * R
- Note: empirically equivalent for typical structures, but semantically different

### 4. Atom selection
- jamaliki: standard protein residues only (32824 atoms for 5A1A)
- MapQ: all non-H atoms (33696 atoms for 5A1A)

### 5. Points per shell
- jamaliki: exactly num_pts (8) per shell, fixed 168-dim correlation
- MapQ: all valid points from successful retry, variable-width masked correlation

### 6. Averaging scope
- MapQ reports Q-scores for protein/nucleic acid atoms, excluding waters
- Including waters inflates the average (waters have high Q ~ 0.79 due to
  spherical density around single atoms)

## Contribution of each difference (EMD-2984/5A1A)

| Cumulative fixes | Q (all) | Q (protein) | Diff from MapQ |
|-----------------|---------|-------------|----------------|
| jamaliki original | 0.644 | - | +4.7% |
| + Saff-Kuijlaars + 0.9R + all atoms | 0.624 | - | +1.4% |
| + All valid points per shell | 0.619 | 0.615 | +0.7% / +0.1% |
| Exclude waters from average | - | 0.615 | +0.1% |
