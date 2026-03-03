# Q-score Implementation Differences: jamaliki/qscore vs MapQ

## Summary

Five implementation differences in jamaliki/qscore explain the systematic +4% offset
from MapQ reference values. When all five are addressed, agreement improves from
+4.7% to +0.4% (mean across test entries).

## Differences

### 1. Epsilon stabilization (critical for low-variance maps)
- **jamaliki**: `eps=1e-6` in Pearson denominator
- **MapQ**: no epsilon (`numpy.corrcoef` / denominator > 0 check)
- **Impact**: catastrophic for maps with std < 0.001 (EMD-72359: 0.08 vs 0.38)
- **Fix**: reduce epsilon to 1e-12 or use denominator > 0 guard

### 2. Sphere sampling method
- **jamaliki**: random (`np.random.randn`), non-deterministic, retries until 8 valid
- **MapQ**: Saff-Kuijlaars spiral, deterministic, cached
- **Impact**: per-atom Q-score variance up to 0.25 between runs

### 3. Neighbor exclusion criterion
- **jamaliki**: point must be the closest to its originating atom (KDTree nearest query)
- **MapQ**: nearest non-self atom must be farther than 0.9 * R
- **Impact**: different subsets of sphere points are selected

### 4. Atom selection
- **jamaliki**: standard protein residue types only, standard atom names only
- **MapQ**: all non-hydrogen atoms including HETATM
- **Impact**: ~2-3% fewer atoms in jamaliki (e.g., 32824 vs 33696 for 5A1A)

### 5. Points per shell in correlation
- **jamaliki**: exactly `num_pts` (8) values per shell, 168-dim correlation
- **MapQ**: ALL valid points from the successful retry attempt, variable-width masked correlation
- **Impact**: shells at larger R with many valid points contribute more information

## Verification

| Entry | jamaliki | MapQ-matching | MapQ ref | Improvement |
|-------|----------|---------------|----------|-------------|
| EMD-2984 (2.2Å) | 0.644 (+4.7%) | 0.619 (+0.7%) | 0.615 | 4.0% closer |
| EMD-72359 (3.9Å) | 0.085 (-77.5%) | 0.376 (+0.1%) | 0.376 | 77.4% closer |

## Recommendation

For IHMValidation integration, either:
1. Apply all five fixes to jamaliki/qscore (substantial refactor)
2. Use jqscore (github.com/aozalevsky/jqscore) which already implements MapQ-matching logic
3. Implement a new pure-Python module incorporating these findings
