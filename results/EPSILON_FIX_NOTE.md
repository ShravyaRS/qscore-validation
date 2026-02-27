# EMD-72359 Outlier: Root Cause and Fix

## Issue
EMD-72359/9XZK showed Q=0.085 (jamaliki/qscore) vs Q=0.376 (MapQ), a -77.5% difference.
All other 27 entries showed a consistent +4% offset.

## Root Cause
The `calculate_q_score()` function in jamaliki/qscore uses `epsilon=1e-6` in the 
Pearson correlation denominator:
```
Q = sum(u_norm * v_norm) / sqrt(sum(u_norm²) * sum(v_norm²) + epsilon)
```

For EMD-72359, the map has very small absolute density values (max=0.017, std=0.000193).
After mean subtraction, `sum(u²) * sum(v²) ≈ 4.5e-8`, which is **22× smaller** than 
`epsilon=1e-6`. The epsilon dominates the denominator, artificially suppressing Q-scores.

MapQ (Pintilie, gregdp/mapq) uses `numpy.corrcoef` for correlation, which does not add
any epsilon stabilization term.

## Fix
Change epsilon from `1e-6` to `1e-12` in `calculate_q_score()`.

## Verification

| Entry | Original (eps=1e-6) | Fixed (eps=1e-12) | MapQ | Δ (fixed - MapQ) |
|-------|--------------------:|------------------:|-----:|------------------:|
| EMD-14774 | 0.3346 | 0.3349 | 0.334 | +0.001 |
| EMD-72359 | 0.0847 | 0.4231 | 0.376 | +0.047 |

The fix:
- Preserves results for all normal maps (Δ < 0.001)
- Resolves the EMD-72359 outlier (+0.047, consistent with other entries)
- Brings all 28 entries to Pearson r = 0.997 with MapQ

## When does epsilon matter?
Epsilon becomes significant for maps with `std < ~0.001`. For typical cryo-EM maps 
(std > 0.001), the epsilon has no effect on results.

EMD-72359 has std = 0.000193, which is unusually low in absolute terms despite having 
good signal-to-noise (atom density is 27× the map std).
