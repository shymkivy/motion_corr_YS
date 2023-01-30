# motion_corr_YS

for movie preprocessing

use "preprocess_mov.m" for processing

used motion correction and bidi shifts fix shifts from suite2p matlab repository as starting point, and tuned to work better with noisy movies and my datasets
  	https://github.com/cortex-lab/Suite2P

bidi changes, for large shifts
	transform the image with interpolation into shifted sinusoidally transformed coordinate system, since res galvo scans sinusoidally (required for large shifts)

Moco changes, to work with noisy movies
	created an iterative process, where first iterations are done with large xyz smoothing, iteratively applied and slowly reduced
	added regularization to both steps of smoothing, to penalize large shifts
	added cap on high signal pixels, to prevent large correlation values from single bright sources
	Nonrigid still not optimized, so only use rigid
