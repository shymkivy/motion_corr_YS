# motion_corr_YS

for movie preprocessing

Use "preprocess_mov.m"

Started with Suite2P-matlab motion correction and bidi shifts fix, because they are really good and fast, and made some adjustments work better with noisy movies
* https://github.com/cortex-lab/Suite2P

Bidi changes, for large shifts
* transform the image with interpolation into sinusoidal coordinate system during shifts, an then back, since res galvo scans sinusoidally (mostly required for large shifts)

Moco changes, to work with noisy movies
* created an iterative process, where first iterations are done with large xyz smoothing, iteratively applied and slowly reduced
* added regularization to both steps of motion correction, to penalize large shifts
* added cap on high signal pixels, to prevent large correlation values from single bright sources
* Nonrigid still not optimized, so only use rigid
