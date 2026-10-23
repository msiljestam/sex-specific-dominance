# SSdomSL

MATLAB simulation of trait and regulator evolution under sex-specific selection, with autosomal or X-linked inheritance.

## Usage

Requires Statistics and Machine Learning Toolbox; tested with MATLAB R2025b. Place the repository files together on the MATLAB path or open the directory in MATLAB.

```matlab
[n_mean, ne_mean, xr_mean] = Main;
```

Simulation parameters can be supplied as name/value pairs, for example `Main('sif', 1, 'sim', 1)`. See the parameter section in `Main.m` for defaults and descriptions. The default run matches panel ii of the manuscript figure labelled `fig:sim`, with `sif = 1`, `sim = 1`, `x0 = 1`, and 200,000 mutation-input generations.

The mutation-free cooldown is always `Tcd = min(2e5,round(T/2))`, where `T` is the number of mutation-input generations. For the default `T = 2e5`, this gives `Tcd = 1e5`.

The fixed runtime and plotting settings `prpl`, `prog`, `t0`, `Nlim_rel`, `pl_res`, `Tcd`, and `fs` are defined immediately after `loadinput` in `Main.m`; to change them, edit that section rather than passing name/value arguments.

Plotting uses `t + t0` to accommodate logarithmic time axes; `t0 = 100` also delays intermediate figure updates until after generation 100, although a final plot is always made. `pl_res` targets 1,000 sampled generations and is capped automatically for short runs.

## Outputs

The outputs are the three elements of `mean(nxv,1)`, averaged over the final `min(Tcd,1000)` generations:

- `n_mean`: mean number of distinct trait values with positive allele-copy counts.
- `ne_mean`: mean effective number of trait values (inverse Simpson diversity).
- `xr_mean`: mean retained trait range, `max(xx) - min(xx)`.

Counts are pooled by identical trait value across females and males.

## Conventions

`Nmax` is census size in individuals. Regulators are stored in log10 units, so affinity is `10.^d`. Independently evolving regulators are bounded by `+/-abs(dr)/2`; with `mud = 0`, they follow `d = x.*dr`. `delta_x` is the expected mutational step size; `delta_d` uses natural-log affinity units before scaling and boundary clipping.

Initialization requires an expected count of at least 100 individuals per sex.

## Files

`Main.m` runs the simulation. `SimPlot.m` assembles the figure and contains the local plotting helpers. `AllelePlot.m` draws the tree plot and dominance inset.

## Licence and citation

Copyright (c) 2026 Mattias Siljestam. Released under the [MIT licence](LICENSE).

Please cite the software version used and the accompanying manuscript, *The coevolution of sex-specific dominance and allelic diversity at sexually antagonistic loci*. Software author details are provided in [CITATION.cff](CITATION.cff).
