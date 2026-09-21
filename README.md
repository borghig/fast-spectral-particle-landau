# Fast spectral particle method for the Landau equation

MATLAB code for the numerical examples in *A fast spectral particle method
for the Landau equation*. The implementation treats the spatially homogeneous,
two-dimensional equation with Maxwell-molecule interactions. The type-II blob method is included as the paper's comparison
baseline.

## Quick start

Open this folder in MATLAB and run:

```matlab
test_spectral_particle_BKW
```

Edit `M`, `N`, and `tmax` in the settings section at the top of
`test_spectral_particle_BKW.m`. Their defaults are 64 modes per axis,
80² particles, and final time 8. Choose M = 32, 48, 64, or 128, a perfect
square for N, and a positive multiple of 0.1 for `tmax`.

The remaining defaults use particles initialized on [-5,5]², Fourier domain
[-10,10]², explicit Euler with time step 0.01, and Gaussian filter variance
0.64 h². Diagnostics are recorded every 0.1 time units on a 120 × 120 grid.
The animation shows the density, particles, and a density cross-section
overlaid with the exact BKW solution at the same time and slice location.
The error panel shows relative L², momentum, and energy errors on a logarithmic
scale. The final relative L² error is printed in the Command Window.

For a custom run:

```matlab
setup_paths;
prm = initialization(48, 80^2, 'BKW');
prm.epsi = 0;
prm.G = 1;
prm.dt = 1e-3;
prm.timeint = 'RK4';
[out, prm] = method_spectral_particle(prm);
out = compute_error(out, prm);
```

The supported initial conditions are `BKW`, `Trubnikov`, and `BumpOnTail`.
The particle count must be a perfect square. Precomputed collision modes
are supplied for M = 32, 48, 64, 128. If changing `tmax` or `dt_error`, update
`prm.times_error = 0:prm.dt_error:prm.tmax` as well; diagnostic times must lie
on the time-stepping grid. The bump-on-tail generator performs its own
density reconstruction because a time-dependent exact density is unavailable.

## Paper figures

Run `setup_paths` once, then call the functions below. Every generator finds
its dependencies and output directory relative to its own location.
Saved data and the ten original numerical figure PDFs are under
`figures generation/results/`.

| Paper figure | Function | Output PDF(s) |
| --- | --- | --- |
| 2: BKW reconstruction at t = 0 | `experiment_t0` | `experiment_t0_fig1.pdf` |
| 3a: BKW error over time; 4: runtime and accuracy | `experiment_MvsN` | `experiment_MvsN.pdf`, `experiment_MvsN_cost_accuracy.pdf` |
| 3b: BKW error at t = 8 | `experiment_MvsN_fixed_time` | `experiment_MvsN_fixed_time.pdf` |
| 5: conservation and entropy | `experiment_conservation_entropy` | `experiment_conservation_entropy.pdf` |
| 6: filtering, M = 48 and 128 | `experiment_filter_comparison` | `experiment_filter_M48_N80.pdf`, `experiment_filter_M128_N80.pdf` |
| 7: Trubnikov relaxation | `experiment_trubnikov_temperature` | `experiment_trubnikov.pdf` |
| 8: bump-on-tail relaxation | `experiment_bump_on_tail` | `experiment_bump_on_tail_density.pdf`, `experiment_bump_on_tail_marginal.pdf` |

The generators use the same sequence: settings, load or compute, then plot
and export. Calling one without arguments reuses saved runs and computes
missing data. Passing `true` recomputes that experiment's simulations, for
example `experiment_bump_on_tail(true)`. The instantaneous `experiment_t0`
always recomputes and takes no argument. Running a generator exports to its
listed PDF filename.

`experiment_blob_only` supplies the BKW baseline required by the comparison
figures; dependent generators call it automatically. To recompute that
baseline, call `experiment_blob_only(true)` separately. Blob runs use explicit
Euler with time step 0.01 and can be substantially slower than the spectral
particle runs.

## Files and dependencies

- `method_spectral_particle.m`: spectral particle time integration.
- `method_blob.m`: type-II blob baseline.
- `test_spectral_particle_BKW.m`: the single demonstration.
- `aux_functions/`: initialization, flux and velocity evaluation, density
  reconstruction, diagnostics, and the required Gaussian kernels.
- `modes/`: precomputed collision-kernel Fourier coefficients.
- `figures generation/`: only the experiments and baseline needed by the paper.
- `finufft/`: the required CPU MATLAB interface and Apple silicon MEX binary,
  with upstream license notices.

The original code targets MATLAB R2026a. The bundled FINUFFT runtime is for
Apple silicon macOS. For another platform, install a compatible
[FINUFFT MATLAB CPU interface](https://finufft.readthedocs.io/en/latest/matlab.html)
and add its interface and MEX folders to the MATLAB path before calling
`setup_paths`. Only type-1 and type-2 two-dimensional transforms are used,
with tolerance 1e-13. No MATLAB toolbox functions are used by the project code.

## Citation

If you use this code, please cite:

Giacomo Borghi and Lorenzo Pareschi. *A fast spectral particle method for the Landau equation*. Manuscript, 2026.

```bibtex
@misc{borghi2026fast,
  title         = {A fast spectral particle method for the {Landau} equation},
  author        = {Borghi, Giacomo and Pareschi, Lorenzo},
  year          = {2026},
  archivePrefix = {arXiv},
  eprint        = {},
  url           = {}
}
```
