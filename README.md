# README

This is a fork of the "WRFlux" available at: [https://github.com/matzegoebel/WRFlux](https://github.com/matzegoebel/WRFlux) which is a fork of the "Weather Research and Forecast model (WRF)" available at: [https://github.com/wrf-model/WRF](https://github.com/wrf-model/WRF).

**WRFstats** allows to output time-averaged quantities, and 2nd and 3rd order products of quantities to be able to compute variances, and covariances in post-processing. Other small modifications have been made such as :
- A bug fix in tslist
- A feature to allow the time series maximum height to be different in each domain
- A feature ('fine\_LU') that initialize soil properties by interpolating from parent domain (in case a nested domain starts later than its parent)
- A feature to chose the default soiltype that is used in land cells where the soil category is not defined.
- A bug fix (probably not the best way) [https://github.com/wrf-model/WRF/issues/1947](https://github.com/wrf-model/WRF/issues/1947)

## Usage

### Branchs

The diverse branches corresponds to the different implementation (stats for statistics), and run is the merge of all branches used for the paper "Observation and modeling of an offshore atmospheric undular bore generated during sea-breeze initiation near a peninsula" (commit 7b8cf696)

### Namelist options

During the model run, fluxes, tendencies, and budget variables are averaged over time.
The online calculations can be controlled in the namelist file or the registry file `Registry/registry.wrfstats`. The calculations do not affect the model evolution.

The namelist variable **`output_stats`** is added. It can take the following values :
- 0 (default): no statistics is computed
- 1: return only 1st order statistics (mean)
- 2: return 1st and 2nd order statistics (mean, non-centered variances, non-centered covariances)
- 3 (or more): return 1st, 2nd, and 3rd order statistics

All variables are output to the auxiliary output stream `auxhist24`. The output interval can be set with the namelist variables `auxhist24_interval_m` and `auxhist24_interval_s`. The averaging period is equal to the output interval

### Time-averaged variables

The averaged variables can be found in `Registry/registry.wrfstats`, and the calculation is found in `dyn_em/module_statistics.F`. To add new variables, just add them in `Registry/registry.wrfstats`, and add their calculation in `dyn_em/module_statistics.F`

### Post-processing

Variances and covariances can be computed in post-processing as 

var(u) = u2_avg - u_avg\*u_avg

covar(u, v) = uv_avg - u_avg\*v_avg

## Installation

This repository contains a complete, standalone version of WRF. Since it is a fork of WRF, the whole history of WRF's master branch is included as well. Thus, you can simply merge your changes with the changes that WRFstats introduces using `git`.

## Implementation

### Compute the average

The calculation of time-average values is basic:

$\overline{x}\_n = \sum_{i=0}^{i=n}(x\_i) / n$

Then :

$\overline{x}\_0 = 0$

$\overline{x}\_{n+1} = (n\*\overline{x}\_n + x\_{n+1}) / (n+1))$

It is computed in `dyn_em/module_statistics.F`. All the variables are averaged in the cell centers as in the timeseries feature. Velocities are interpolated before being averaged. This is mandatory to compute covariances. 

The calculation of the water vapor density $\rho_d$ is based on <cite>Skamarock et al. (2019)</cite>, Equation 2.16. The only averaged mixing ratio is the water vapor. The TKE\_AVG variable is the time-averaged of either the TKE\_PBL (from MYJ PBL model for example) or the TKE variable (subgrid TKE from TKE1.5 model). For some variables (e.g. pressure, temperature), the perturbation is averaged rather than the total quantity 

### List of modified files

The following WRF source code files have been modified:
- README\_WRFlux.md (just changed the name)
- Registry/Registry.EM_COMMON (add TKE to the history)
- Registry/registry.em_shared_collection
- Registry/registry.wrfstats (new)
- dyn_em/Makefile (add dependency with module_statistics)
- dyn_em/depend.dyn_em (add dependency with module_statistics)
- dyn_em/module_initialize_real.F
- dyn_em/module_statistics.F (new)
- dyn_em/solve_em.F (stats)
- dyn_em/start_em.F (stats)
- main/depend.common (add dependency with module_statistics)
- share/interp_fcn.F (fine_lu)
- share/mediation_integrate.F
- share/module_check_a_mundo.F
- share/start_domain.F (fine_lu)
- share/wrf_timeseries.F (fix bug)
- wrftladj/depend.wrftladj
- wrftladj/solve_em_ad.F (stats)
- wrftladj/solve_em_tl.F (stats)

## Caveats and limitations

* WRFlux can increase the runtime by about 25 %. No test have been realized for WRFstats but the writing time might be highly impacted
* The tool does not work with adaptive timestepping

## Tests

The online calculation have been compared to time series output


## Contributing

Feel free to report [issues]() on Github.
You are also invited to fix bugs and improve the code yourself. Your changes can be integrated with a [pull request]().

## Acknowledgments

Thanks to Matthias Göbel who provided the code used as a starting point for WRFstats, and WRF developers and support for their amazing work.

