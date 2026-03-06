# Direct Verification from MapQ Source Code (gregdp/mapq)

All key implementation details confirmed by reading gregdp/mapq/qscores.py:

## 1. Saff-Kuijlaars Spiral (SpherePts, line 4760)
```python
def SpherePts ( ctr, rad, N ) :
    for k in range ( 1, N+1 ) :
        h = -1.0 + ( 2.0*float(k-1)/float(N-1) )
        phis.append ( acos(h) )
        thetas.append ( 0 if k == 1 or k == N else
                        (thetas[k-2] + 3.6/sqrt(N*(1.0-h**2.0))) % (2*pi) )
```
CONFIRMED: Deterministic Saff-Kuijlaars spiral, identical formula.

## 2. 0.9R Exclusion (line 698, 901)
```python
outRad = RAD*0.9
# comment: "the 0.9 factor is a small adjustment so that the atom in
# question is not found as a 'nearby' atom"
```
CONFIRMED: 0.9 * R threshold for neighbor exclusion.

## 3. Retry Loop (line 712)
```python
for i in range (0, 50) :
    outPts = SpherePts ( at.coord(), RAD, npts+i*2 )
    ...
    if at_pts_i >= npts or show :
        pts.extend ( at_pts[0:at_pts_i] )  # ALL valid points kept
        break
```
CONFIRMED: 50 retries, increasing by 2 each time, ALL valid points from
successful attempt are kept (not just first num_pts).

## 4. No Epsilon in Correlation (line 1541)
```python
olap, CC, CCm = FitMap.overlap_and_correlation ( d_vals, g_vals )
qscore = CCm
```
CONFIRMED: Uses Chimera's FitMap.overlap_and_correlation which computes
correlation about the mean (CCm) without any epsilon stabilization term.

## 5. Reference Gaussian (line 1504-1506)
```python
A,B = maxD - minD, minD
gv = A * numpy.exp ( -0.5 * numpy.power(RAD/sigma,2) ) + B
g_vals = numpy.append ( g_vals, numpy.ones([len(pts),1]) * gv )
```
CONFIRMED: Same formula as jamaliki. Height = maxD-minD, offset = minD.
One g_val entry per valid point (not averaged per shell).

## 6. Flat Vector Correlation (line 1541)
The correlation is computed on the full flat d_vals and g_vals vectors,
which contain ALL valid points from ALL shells concatenated. This is NOT
a per-shell-averaged 21-dim correlation.
