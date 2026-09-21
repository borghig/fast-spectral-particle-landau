# FINUFFT runtime

This folder contains the original CPU MATLAB wrappers needed by the spectral
particle method: `finufft2d1`, `finufft2d2`, their shared plan and validation
helpers, and the accompanying Apple silicon `finufft.mexmaca64` binary.

The files are copied without modification from the original bundled FINUFFT
checkout. Unused interfaces, examples, tests, build files, and source trees
are excluded from this submission. `LICENSE`, `NOTICE`, and the bundled
dependency license files are retained.

For another platform or to rebuild, use the
[upstream FINUFFT project](https://github.com/flatironinstitute/finufft) and its
[MATLAB installation instructions](https://finufft.readthedocs.io/en/latest/matlab.html).
Add the matching MATLAB interface and compiled MEX directory to the path.
