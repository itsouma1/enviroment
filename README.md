# climate

A dependency-free C++17 implementation of the retrieval formulas, time-series
decomposition and trend tests described in the NASA hackathon research notes:
sea surface temperature (NLSST), atmospheric CO2 curve fitting, NDVI, NASA Team
sea ice concentration and extent, GPI/GPCP precipitation merging, STL
decomposition, Mann-Kendall, Sen's slope, and the deterministic
intent -> catalog -> analysis pipeline.

Only the C++ standard library is used. Nothing is linked, and no numerical
library, JSON library or CSV library is required.

## Build

```sh
make            # builds ./climate and tools/generate_demo
make test       # builds and runs the self-checking test program
make data       # writes the synthetic demo datasets into data/
make demo       # query "is the arctic changing?" against data/
```

Requirements: a C++17 compiler (g++ or clang++) and `make`.

## Layout

| Path | Contents |
|---|---|
| `src/stats.*` | descriptive statistics, normal CDF/quantile, autocorrelation |
| `src/linalg.*` | Gaussian elimination and Householder-QR least squares |
| `src/fft.*` | radix-2 and Bluestein FFT, Thoning-style FFT low-pass filter |
| `src/loess.*`, `src/stl.*` | local regression and the full STL inner/outer loop |
| `src/harmonic.*` | polynomial + harmonic regression, derivative, growth rate |
| `src/trend.*` | Mann-Kendall (tie corrected, Hamed-Rao), Sen's slope with CI, TFPW, seasonal MK |
| `src/series.*` | gap filling, phase climatology, standardised anomalies |
| `src/retrieval.*` | NLSST, NASA Team, NDVI, GPI, GPCP merge, energy/carbon/AMOC diagnostics |
| `src/catalog.*` | the variable catalog with process and region tags |
| `src/intent.*` | deterministic intent parser, JSON IR ingestion, whitelisted compiler |
| `src/analysis.*` | the seasonality -> trend pipeline for ordinary and CO2 series |
| `src/report.*` | text, Markdown and JSON reports |
| `src/main.cpp` | command-line interface |
| `tests/test_all.cpp` | formula and pipeline checks |
| `tools/generate_demo.cpp` | generator for the synthetic datasets |

## Commands

```sh
./climate catalog
./climate catalog --process sea_ice --format json

./climate query "is the arctic changing?"
./climate query "how has sahel rainfall changed since 1980?" --data data
./climate intent --file ir.json

./climate analyze --file data/sst_monthly.csv --variable sst --units degC
./climate analyze --file data/sst_monthly.csv --format json --out sst.json
./climate co2 --file data/co2_monthly.csv

./climate retrieve sst --file data/sst_matchups.csv --out data/nlsst_model.csv
./climate retrieve sst --coeffs data/nlsst_model.csv --bt11 300.2 --bt12 299.3 \
        --first-guess 300.0 --zenith 35
./climate retrieve ndvi --file data/reflectance_pairs.csv --red red --nir nir
./climate retrieve seaice --file data/sea_ice_grid_2024_09.csv
./climate retrieve seaice-pixel --tie-points data/nasa_team_tiepoints.csv \
        --tb19v 240 --tb19h 228 --tb37v 242
./climate retrieve precip --file data/precip_estimates.csv --cold-fraction 0.35 \
        --water-vapour 4.2 --a1 0.1 --a2 3.0

./climate calibrate nlsst --file data/sst_matchups.csv --out data/nlsst_model.csv
./climate calibrate seaice --tie-points data/nasa_team_tiepoints.csv --out data/nasa_team_model.csv

./climate physics toa --forcing 3.7 --feedback -1.2 --delta-t 1.2
./climate physics surface --albedo 0.2 --shortwave 200 --emissivity-air 0.9 \
        --longwave 300 --emissivity-surface 0.98 --temperature 288 \
        --latent 80 --sensible 40 --ground 5
./climate physics carbon --start-ppm 400 --end-ppm 410 --start-year 2000 --end-year 2010
./climate physics amoc --file data/amoc_velocity_26n.csv --depth 500
./climate physics feedback --planck -3.2 --water-vapour 1.8 --lapse-rate -0.6 \
        --albedo 0.35 --cloud 0.4
./climate physics forcing-co2 --ppm 420 --reference-ppm 280
./climate physics mie --size-parameter 5 --real 1.33 --imag 0.01
./climate physics aod --number-concentration 1e8 --radius-um 0.1 --path-length 1000
./climate physics wind-stress --u10 10 --v10 0
./climate physics ekman --tau-x 0.1 --tau-y 0 --latitude 30
./climate physics geostrophic --density 1025 --dp-dx 1e-4 --dp-dy 0 --latitude 30
./climate physics thermal-wind --density 1025 --drho-dx 0 --drho-dy 1e-6 --latitude 30

./climate earth two-layer --years 150 --forcing 3.7
./climate earth two-stream --tau 5 --omega 1 --asymmetry 0.85 --zenith 60
./climate earth equatorial-wave --k 2e-7 --c 30 --mode 0
./climate land penman-monteith --net-radiation 100 --delta 1 --gamma 0.066
./climate land farquhar --vcmax 100 --ci 30 --gamma-star 4 --kc 40 --oi 21 --ko 21000 --j 150
./climate spatial getis --values "1,2,3" --weights "1,1,1;1,1,1;1,1,1" --index 1

./climate atmos cape --file data/parcel_profile.csv
./climate atmos photolysis --file data/photolysis_bands.csv
./climate atmos falloff --temperature 300 --density 1e19
./climate ocean mixed-layer --file data/n2_profile.csv
./climate ocean osborn --gamma 0.2 --epsilon 1e-8 --n2 1e-4
./climate ocean waves --wind 10 --omega 1
./climate ocean hibler --p-star 2e4 --thickness 2 --concentration 0.9 --c 20

./climate retrieve index --file data/reflectance_pairs.csv --red red --nir nir \
        --index evi --blue blue

./climate data list
./climate data inspect --file FILE.nc
./climate data export --file FILE.nc --variable sst --out sst.csv --mode timeseries
./climate data search --dataset sea_ice_extent --start 2024-01-01 --end 2024-12-31
./climate data fetch --dataset sea_ice_extent --start 2024-01-01 --end 2024-12-31 --max 12
./climate data convert --dataset sst --out sst.csv granule1.nc granule2.nc
```

Run `./climate help` for the complete option list.

## Implemented formulas

**SST (NLSST).** `SST = a0 + a1 BT11 + a2 (BT11-BT12) T_sfc
+ a3 (sec(theta)-1)(BT11-BT12) + a4 mirror + a5 theta* + a6 theta^2`.
Coefficients are estimated by Householder-QR least squares from
satellite-buoy matchups (`calibrate nlsst` / `retrieve sst --file`), or loaded
from a coefficient CSV. The retrieval is applied per footprint with
`retrieve sst --bt11 ...`.

**CO2 (NOAA/CMDL).** Polynomial of degree *k* (default 2) plus *n* annual
harmonics, fitted with iterative rejection of points beyond 3 residual standard
deviations, exactly as in Thoning, Tans & Komhyr (1989). Harmonic amplitudes are
`sqrt(A_n^2 + B_n^2)` and phases `atan2(B_n, A_n)`. The growth rate is the
derivative of the fitted trend (polynomial part only, so the seasonal cycle does
not contaminate it) and is also converted with `k = 2.124 PgC/ppm`. An
FFT low-pass filter with configurable cutoff is available for the
short-term/long-term separation.

**NDVI.** `(rho_NIR - rho_RED) / (rho_NIR + rho_RED)`, with an undefined result
(a non-finite value) when both reflectances are zero.

**Sea ice.** The NASA Team algorithm is implemented in two forms:

* direct inversion of the physical three-constituent mixture model for the
  open-water, first-year and multiyear fractions using a 3-channel tie-point
  matrix plus the closure `C_OW + C_FY + C_MY = 1`;
* the empirical polarization/gradient form
  `C = (a0 + a1 PR + a2 GR + a3 PR GR) / (1 + d1 PR + d2 GR + d3 PR GR)`,
  with PR(19) and GR(37V/19V) and coefficients calibrated by least squares over
  the mixture simplex generated from the tie points.

Extent uses the 15% threshold: `SIE = sum A_cell * 1[SIC >= 0.15]`, area uses
`SIA = sum SIC * A_cell` with spherical cell areas
`R^2 dlambda (sin(lat2) - sin(lat1))`.

**Precipitation.** GPI = 3 mm/h times the fraction of IR pixels colder than
235 K, optionally scaled by `exp(a1 (w - a2))`, then a GPCP inverse-error-variance
merge with weights `1/sigma_i^2`.

**Coupled-system diagnostics.** `N = F + lambda dT`,
`Rn = (1-alpha) Rs + eps_a Rl - sigma eps_s Ts^4` with the
`Rn = LE + H + G` residual, the carbon mass-balance growth rate, and the
zonally/depth-integrated AMOC streamfunction `Q(y, z)` in Sverdrup.

**Radiative transfer and dynamics.** Feedback decomposition
`lambda = lambda_0 + sum_i lambda_i` with equilibrium warming `-F/lambda` and
ECS `-F_2x/lambda`; logarithmic CO2 forcing `5.35 ln(C/C0)`; Mie scattering for
a homogeneous sphere (downward `D_n` recurrence, Riccati-Bessel functions, Qext,
Qsca, asymmetry `g`, single-scattering albedo); aerosol optical depth
`N pi r^2 Qext L`; bulk wind stress `rho_a Cd |U| U`; Ekman transport
`(tau_y, -tau_x)/(rho f)` and pumping `curl(tau)/(rho f)`; geostrophic balance
`f u = -(1/rho) dp/dy`, `f v = (1/rho) dp/dx`; thermal wind
`f du/dz = (g/rho) d rho/dy`, `f dv/dz = -(g/rho) d rho/dx`; and the SAVI and
EVI vegetation indices.

**Advanced climate sensitivity.** Two-layer energy balance
`C dT/dt = F + lambda T - eps gamma (T - Td)`, `Cd dTd/dt = gamma (T - Td)`,
integrated with classical RK4; ocean heat uptake `gamma (T - Td)`, its efficiency
`OHU / dT`, transient response `F_2x / (-lambda + eps gamma)`, and the
pattern-corrected sensitivity `ECS = -F_2x / (alpha + alpha')`. Carbon-cycle
feedback decomposition `dC = beta dC_atmos + gamma dT` for land and ocean, the
airborne fraction `dC_atmos * 2.124 / emissions`, and permafrost decay
`dC/dt = -k(T) C` with `k(T) = k_ref Q10^((T-T_ref)/10)`.

**Atmospheric radiation and turbulence.** Planck radiance
`2hc^2/lambda^5 / (exp(hc/lambda kT) - 1)` and its inverse for brightness
temperature; Schwarzschild optical depth, transmittance, weighting function
`(1/g) k q exp(-tau)` and contribution function; a two-stream Eddington solver
with direct beam, surface albedo and finite-difference solution of the coupled
flux equations (energy closure checked by `make test`); shallow-water and
equatorial wave dispersion `omega^3/c^2 - omega(k^2 + (2n+1)beta/c) - beta k = 0`;
Kolmogorov `phi_n = 0.033 Cn^2 k^(-11/3)`, dissipation rate, microscale,
`E(k)` spectrum and the stable boundary-layer depth ratio.

**Land surface and biosphere.** Net radiation, bulk sensible heat, aerodynamic
resistance, Penman-Monteith latent and sensible heat (which sum to `Rn - G` by
construction), soil-moisture surface resistance, Farquhar Rubisco- and
electron-transport-limited assimilation, Ball-Berry, Medlyn and Leuning
stomatal conductance, LPJ NPP and maintenance respiration, LAI from leaf carbon,
xylem cavitation conductance, the soil water balance and soil thermal inertia.

**Spatial statistics.** Getis-Ord Gi*, local Moran's I, the Geo-Detector
q-statistic, the Walker field-significance threshold, Benjamini-Hochberg
q-values and rejections, and the Kulldorff Poisson scan log-likelihood ratio.

**Atmospheric chemistry and aerosols.** CAPE and CIN by trapezoidal integration of
`g (Tv_p - Tv_e)/Tv_e`; photolysis `J = int phi sigma F dlambda`; the
low-pressure-limit falloff constant with the `0.6^Z` broadening factor;
aqueous sulfate oxidation; the free-molecular coagulation kernel and
coagulation sink; the ion balance steady state (and its
`sqrt(q/alpha)` limit); and the ACI susceptibility decomposition
`ERFaci = (dCF/dlnAI * SW + CF * dSW/dlnAI) dlnAI`.

**Ocean physics.** Buoyancy frequency, the `N^2`-weighted mixed-layer depth,
Ekman velocity and barotropic compensation, PV and GM eddy diffusivities, the
internal-wave dispersion `omega^2 = f^2 + N^2 k^2/m^2`, gradient Richardson
number, Osborn diffusivity, shear dissipation, the NORi piecewise viscosity,
Martin export, air-sea CO2 flux, DIC, bulk latent/sensible heat fluxes, iron
co-limitation and ligand speciation, mCDR potential, the Pierson-Moskowitz
spectrum, Doppler-shifted frequency and Hibler ice pressure.

**Closure upgrades.** Convective updraft budget with the diagnosed drag length
`H = r_d sqrt(a_u)` and its `d/dH` sensitivity; mass-based free-molecular
coagulation and the log-normal moment decay; IDEMIX `mu0`, quadratic wave
dissipation and its `d/dE`; the KPP bulk Richardson number with `dRib/dVt`;
BBM damage-dependent viscosity with its derivative, Mohr-Coulomb yield, damaged
elasticity and the three friction regimes; and the martingale-pressure
acceleration under the strong hydrostatic hypothesis. TEOS-10 is provided as
generic Gibbs evaluators (pure-water 8x7 coefficients and saline terms) plus a
sound-speed routine that differentiates any supplied Gibbs function.

**TEOS-10.** The 41-term IAPWS-2009 pure-water table is embedded as an 8 x 7
coefficient array, and `teos10_pure_water_derivatives` returns the full tuple
`(g, g_T, g_P, g_TT, g_TP, g_PP)`. `sound_speed_from_derivatives` uses
`c = g_P sqrt(g_TT / (g_TP^2 - g_TT g_PP))`; at 10 C and 1000 dbar it reproduces
the official GSW value 1463.5311 m/s. The saline evaluator takes the same term
structure but still needs its 64 coefficients; GSW-C embeds them as expanded
polynomial code rather than a data array, so they are not copied blindly.
`climate ocean teos-pure --t 10 --p 1e7` prints the tuple and sound speed.

**Cryosphere (notes 96-112 KB).** `src/cryo.*` implements Glen's rate factor
and flow law, the effective strain/stress invariants, ice viscosity, the SIA
velocity and Weertman-Budd/Coulomb sliding, bulk air/sea-ice stress with the
rotation matrix, the full sea-ice momentum tendency, surface melt energy, the
Stefan front velocity, Stefan number and active-layer depth, three-equation and
simplified basal melt, the plume melt law, subglacial discharge and hydraulic
potential, geothermal/frictional bed melt, calving flux, the phase-field notch
depth, Arthern firn densification with the 550 kg/m^3 rate switch, overburden
pressure, the Herron-Langway gradient, the Nye vertical-velocity profile, the
tipping-point equilibrium and fold thickness, and the scalar Kalman update.
**Advanced cryosphere (notes 112-128 KB).** Damage mechanics is covered by
`damaged_viscosity`, `damaged_glen_strain` with the `(1-D)^-n` enhancement,
`biot_effective_stress`, and `maxwell_elasto_brittle_stress`; the Lagrangian
granular branch by `floe_acceleration`; channelised drainage by
`roethlisberger_channel_rate` and `channel_energy_melt`; phase-field fracture
by `phase_field_energy_density` and `phase_field_evolution` (using the
variational derivative with a supplied Laplacian); plume dynamics by
`plume_flux_derivative` and `plume_momentum_derivative`; probabilistic dating
by `probabilistic_age_increment`. The land-surface formulas in the same window
(net radiation, Monin-Obukhov resistance, Penman-Monteith, soil resistance,
Farquhar, Ball-Berry-Leuning, Medlyn, Leuning, LPJ NPP and maintenance
respiration) were already implemented in `land`.

**Supplied formula gaps.** The Hayhurst invariant is implemented as
`alpha sigma_I + beta sigma_vM + (1-alpha-beta) tr(sigma)` with its invariant
sensitivities and the threshold damage rate `(sigma_H/sigma_th)^r` plus its
derivative. DEM contacts now use the Hertz-Mindlin law: normal stiffness and
viscous damping from contact area, equivalent mass and packing parameter,
Mindlin shear stiffness `K_n/(2(1+nu))`, and the Coulomb cap `mu|F_n|`.
Firn strain softening is the full strain-rate tensor norm divided by the
classical vertical strain, applied as a scale factor to the densification
rate. IceChrono's quadratic cost is evaluated by Cholesky whitening of the
correlation matrix, and its gradient returns `2 G^T P^-1 r` via forward and
back substitution.

**Land and biosphere (notes 128-144 KB).** Added leaf turnover with the
specific-leaf-area longevity law, layered root uptake, the moisture-limited
transpiration power law, the hydraulic-penalty stomatal objective, the
trait-based profit function, growth-temperature acclimation, mesophyll
limitation, CNP colimitation, stoichiometric nutrient demand, daily nitrogen
fixation, canopy air-space exchange, the multilayer scalar tendency, the
exponential soil-production function and its runoff/temperature-dependent
reference rate.

**Ozone/water-vapour feedback and attribution (notes 192-208 KB).** The
Chapman+2 linear damping rates `kappa_O` and `kappa_O3`, the odd-oxygen
tendency, the effective `kappa_Ox`, the SWV feedback `lambda = R/Delta T_s`,
the tropical adjustment `R_atmos - R_stra`, the Jacobian temperature adjustment
`Delta T = -J^-1 Delta H` (Gaussian elimination), the lagged methane-water
trend `2 alpha dCH4/dt`, the cold-trap sensitivity `0.45 ppmv/K`, kernel SARF
summation, and ridge regression `(X^T X + lambda I)^-1 X^T y` with a linear
predictor for the mloz/MLR emulators.

**Radiative attribution protocols (`attrib.hpp`).** One-sided and two-sided
PRP feedback `(R_pert - R_ctrl)/(Ts_pert - Ts_ctrl)`, the fSST ERF as the
perturbed-minus-control mean imbalance, `ERF = IRF + sum(adjustments)`, the
published cloud-adjustment bracket divided by an explicitly supplied
reference-state term, the eight-term surface-temperature linearisation
(six-term core plus ground heat and residual), and the CMIP6 SWV validation
targets (0.25 +/- 0.06 and 0.020 +/- 0.007 W m^-2 K^-1).

**Selectivity, Bayesian structure mapping and edge detection.** The regional
module adds the percentage sign statistic, pooling standard-error ratio,
point-wise Bayesian posterior/log-odds/prior sensitivity, probability-weighted
structure centroids and their sensitivity, the one-sided KS selectivity
statistic and optical-class trend precision. The flow module adds Sobel
gradients, calibrated total edge signal, non-maximum suppression, hysteresis
thresholding, Cohen's d and 3-D connected components. The ocean module adds the
dynamic optical-class mean, static-region mixture and class-switch shift.

**Detectability, quantile clustering and flow geometry.** `power.hpp` gives the
AR(1) variance inflation, MDT at explicit power, detection time, power function,
SNR, sensitivity derivatives and the field-significance inflation. `quantile.hpp`
provides check loss, quantile regression, ALD density and subgradients, Potts
conditional probabilities, Gibbs sweeps and Swendsen-Wang updates. `ensemble.hpp`
builds C0^{d} correlation matrices, Cholesky factors, trend-free AR(1) fields,
sign-test null ensembles and Monte Carlo power. `flow.hpp` implements FTLE and
FTLE fields, the JetLag arc-length descriptor and ridge axis, occurrence and
regime-shift statistics, cyclic zonal smoothing, TempestExtremes-style stitching
and isentropic trend correction. `ocean2.hpp` adds depth warming gradients,
Hobday MHW events, VICA, co-occurrence and independence, subsurface
classification, GMM/PCM posteriors and gradients, K-means and cluster heat
content, REOF-style ridge fingerprints and tracer decomposition.

**Regional sign coherence and Lagrangian objects.** `regional.hpp` implements
the Huth-Dubrovsky sign-difference statistic with binomial and empirical nulls,
effective sample size, critical value and its autocorrelation derivative, and a
normal quantile. It provides the Cai zonal-mean/zonally-asymmetric decomposition,
regional means and sign agreement, GEV location and trend scores, hierarchical
partial-pooling shrinkage and pooled variances, symmetric eigendecomposition,
the North test, varimax rotation, REOF RPC trends and back-projection, the
Weatherhead standard error, detectable trend and detection time with
derivatives, vertical trend-pattern correlations, coherence height, segmented
reversals, and EDW slopes, leverage, latitude correction and lapse-rate trends.
`lagrangian.hpp` implements Hessian eigenvalues, the SCAFET shape index and its
derivatives, grid-aware Gaussian smoothing, object frequency and mean
intensity, nearest-centroid linking and population/along-track trends.

**Variable-specific validation (notes 480-520 KB).** `validation.hpp`
implements the standardised-difference truth test, the quadruple-collocation
closure, LIF, the quadratic sea-ice SIC bias and derivative,
inverse-variance precipitation merging, the hit/missed/false bias sum, the
OCO retrieval+representation error covariance and flux error decomposition,
the ocean-colour `J^T C J` propagation, data/model measurement uncertainty
and the partially-linear CATE estimator.

**Homogenisation and causal attribution (notes 448-480 KB).** `homogen.hpp`
implements the benchmark metrics (absolute bias, FOEX, relative RMSE), the
bivariate ACMANT objective, ensemble homogenisation mean/variance, the
BYS-CUSUM residual statistic and Student-t likelihood, the Singh-Maddala
density, the sparse-jump objective, the cautious selection-bias correction,
SUHI and excursion probability, and the variational lower bound.
`intercal.hpp` adds STAR V5.0 linear/Z terms and TLT synthesis, ACI
susceptibility `beta1 alpha1`, cloud-fraction N_d bias, LWP error variance,
the DML residual estimator, G-computation, natural-effect addition, SNO
scene variance, CLARREO bias variance, TRUTHS total uncertainty, the
random-walk Kalman update, and CycleGAN adversarial/cycle losses.

**SPDE, GTWR and sensor intercalibration (notes 384-400 KB).**
`geostat.hpp` implements the SPDE Matern variance and precision
`tau^2(kappa^2 C + G)^alpha`, GMRF conditional mean/covariance,
squared-exponential GP posterior mean and variance, the scaled MAD envelope
statistic and p-value, the GTWR space-time distance/weights/local estimator
with AICc, the spherical/exponential/Gaussian/Matern (nu=1/2,3/2,5/2)
variograms, and environment-weighted PCA projection. `intercal.hpp`
implements GSICS correction and covariance-aware uncertainty, the MSU/AMSU
nonlinear calibration and drift scaling, double differencing, diurnal and
seasonal harmonics, CWV/SST logarithmic retrievals with state-dependent
sensitivities, SSA hybrid recalibration, CLARREO inverse-variance weights
and upper-bound inflation, the TRUTHS OBCS measurement function, angular
standard deviation, bias-regression slope, and the SNO three-term variance
budget with its t-statistic.

**Metrological uncertainty and intercalibration (notes 352-384 KB).**
`metrology.hpp` implements the covariance-aware LPU in matrix and scalar
forms, the estimator covariance `G Q_y G^T`, the structured decomposition
`A_bar Q_xE A_bar^T + Q_Y + Q_E`, GLS and OLS estimators, Monte Carlo sample
covariance, NDVI partials and correlated uncertainty, simple kriging and its
error covariance `C_cross^T (C_obs+E)^-1 C_cross`, the
`0.5/sqrt(N)` homogenisation and `sigma/sqrt(n)` sampling budgets, bootstrap
mean error, the SNO bias, and the diurnal harmonic reconstruction.

**Advanced magnitude and validation (notes 320-352 KB).** `change.hpp` adds
the piecewise-linear value and rate shift, the non-stationary effective
return level `mu + sigma/xi ([ -ln(1-p) ]^-xi - 1)`, the return-period
derivative `Delta tau/(1-F)^2 dF/dt`, the AR(2) stability sum, the
amplitude/circulation/covariance OHT terms, Fisher information `1/sigma^2`
and information speed, clipped standardisation, local contrast, Max-CUSUM,
higher criticism, the Gaussian covariance log-likelihood ratio (Cholesky
log-determinants), the 1-D optimal-transport cost, and the functional slope
function with its integrated squared norm. `uncertainty.hpp` implements
independent and correlated uncertainty propagation, correlated mean
variance, the retrieval/sampling two-component budget, reduced chi-square,
RMSE, bias and normalised mean bias.

**Spatial detection (notes 288-320 KB).** `spatial.hpp` now adds LISA regime
classification, the Gaussian kernel, GWR local coefficients
`(X^T W X)^-1 X^T W y`, Matern 3/2 covariance, probit changepoint category
probabilities, the Hotelling functional scan, the CRISP objective, LGCP
log-intensity and DIC. Getis-Ord, local Moran, Geo-Detector, Walker, BH and
the Poisson scan were already present. `coherence.hpp` implements the Morlet
CWT and box-smoothed wavelet coherence. `magnitude.hpp` implements quadratic
level/rate/acceleration at the sample end, the linearithmic value, rate and
`beta/t` acceleration, step magnitude, trapezoidal integrated change, return
level shifts and the Newey-West Bartlett long-run variance.

**2026 time-varying inference (notes 256-288 KB).** `regime.hpp` implements
the symmetric alpha-stable increment, the feedback manifold `h(T) = lambda(T)
+ kappa(T)/T` and `dT/dt = F - lambda T - kappa`, logistic and exponential
smooth transitions with the multiple-LSTR mixture, stochastic-volatility
observation/state steps and the `K = -ln(phi)/dt`, `gamma0 = theta/(1-phi)`
estimates, GJR-GARCH, Beta-weighted GARCH-MIDAS, the asymmetric-Laplace
quantile density, seasonal-trend averaging with its covariance standard
error, and BEAST harmonic/step components. `detect.hpp` adds the GEV
log-likelihood, Marchenko-Pastur signal counting, eigenvalue reconstruction
and total least squares; `causal.hpp` provides the panel SEM interaction
design and two-way demeaning. The printed LSTR2 transition `1 + exp(...)`
is outside [0,1], so the bounded `1 - exp(-gamma (x - tau)^2)` is used, and
Beta weights use `k/(K+1)` so endpoint exponents stay finite.

**Causal forests, copulas and HMMs.** `ml.hpp` also provides the honest
CT-H criterion `(1/N_tr) sum tau^2 - (1/N_tr + 1/N_est) sum S^2`.
`copula.hpp` implements the five non-stationary parameterisations (Gumbel,
Joe, Clayton, Huesler-Reiss, Galambos) and their derivatives. `hmm.hpp`
implements scaled Gaussian-emission forward-backward (alpha, beta, gamma,
xi, log-likelihood) and Baum-Welch re-estimation of pi, A, mu and sigma^2.
`climateidx.hpp` exposes `hovmoller_average` and `phase_space_derivative`
as plain arrays for the plotting layer.

**Learners (`ml.hpp`).** A CART regression tree with variance-reduction
splits and impurity importance, a bootstrap random forest with feature
subsampling, the deweathered prediction from a fixed meteorological state,
the honest causal-forest split criterion, and the CATE estimand. These are
the algorithms behind the deweathered-ML and causal-forest diagnostics;
their trained coefficients and splits are data-derived and are not
hard-coded.

**Change detection and extremes (notes 224-240 KB).** GEV and generalized
Pareto CDFs and return levels (parameters may be functions of a covariate),
the Ornstein-Uhlenbeck variance `sigma^2/(2 lambda)` and AR1
`exp(-lambda dt)`, the sliding-window AR1 estimator, the rate-mismatch score,
RSI and CUSUM, optimal fingerprinting `(X^T C_N^-1 X)^-1 X^T C_N^-1 y`, and
the Marchenko-Pastur bounds and density. Mann-Kendall and Sen's slope live in
`trend.hpp`. Climate indices, OHC/EEI, the sea-level budget, quadratic and
broken-linear sea level, the local-linear DLM, groundwater storage, recharge
convolution and depletion, Heston variance and Markov regimes live in
`climateidx.hpp`.

**Land-use change (notes 208-224 KB).** Fraction-weighted albedo mixing, the
annual mean albedo change, kernel and clear-sky `-0.7 Delta alpha SWDW` forcing,
the global conversion sum, `r_c = r_xi/(0.5 LAI)`, canopy Penman-Monteith, the
eight-term surface-temperature decomposition, bookkeeping and stock-difference
carbon accounting, IEF-20, gain-loss stock changes, the four pool tendencies
and the net land-use flux.

**Aerosol optics and forcing (notes 160-192 KB).** A full Mie solver for
homogeneous spheres: logarithmic-derivative downward recurrence, Miller's
downward Riccati-Bessel recursion and the stable chi recursion, returning
`Q_ext`, `Q_sca` and `Q_abs`. It matches the reference `miepython`
implementation to 1e-9 at size parameters 1, 5 and 10 with absorbing and
non-absorbing indices. Around it: cross sections `pi r^2 Q`, column AOD,
single-scattering albedo, Angstrom exponent and spectral scaling, the
Charlson-Sokolik-Toon thin-layer forcing, the two-stream flux derivatives,
polarized reflectance from Stokes Q/U, Kohler activation `N_a w/(w + d N_a)`,
the ACI indices, the forcing decomposition and residual Twomey estimate, and
the chained machine-learning forcing.

**Ozone and stratospheric water vapour (notes 176-192 KB).** Chapman ozone
production and odd-oxygen loss, the NOx catalytic loss rate, the Chapman
equilibrium ratio, methane fractional release, the AoA cubic, the
`H2O + 2CH4` conservation law, methane-driven SWV change, and the column mass
perturbation.

**Dimensional corrections.** Two implemented laws were dimensionally wrong as
supplied and are now fixed rather than copied literally. The Hertz-Mindlin
damping coefficient was given as `2 alpha A xi r_c sqrt(m K_n)/d`, which carries
units of area; it is now the standard `C_n = 2 xi sqrt(m_eq K_n)`, with the
geometric factor retained only in `K_n`. The age-dependent dissolution term
`K tau^sigma` had units of `time^(sigma-1)`; particle age is now normalised by
a reference age, `K (tau/tau_ref)^sigma`, so the factor is dimensionless.

**Five supplied closures.** SPLASH Green-Ampt infiltration
`f = Ks(1 + S(theta_f - theta_i)/F)` with its `d f/dF` sensitivity and the
ponding time `S Ks dtheta / (P(P-Ks))`. Perirhizal hydraulics with the geometry
factor `B(p)`, the effective conductivity between bulk and root potentials, the
conductance `2 pi l_root B k_eff`, and the root-system series conductance with
its sensitivity to `K_prhiz`. The Hybrid-JSBACH4 network is a generic
feed-forward evaluator returning the three outputs and an exact Jacobian by
forward-mode propagation (the published weights are site-specific and were not
released, so callers supply them). Hybrid 4D-Var-EnKF analysis and gain with the
`d/dgamma` sensitivity, and a deterministic EAKF using Cholesky square roots of
the prior and posterior covariances. Multilayer canopy flux, layer assimilation
with the parallel cuticular pathway, and its stomatal sensitivity.

**Critical zone (notes 144-160 KB).** Mineral dissolution constant, silicate
weathering flux by column integration, carbonate weathering with the
runoff-dependent soil pCO2, stream-power erosion, petrogenic and sulfide fluxes,
Michaelis-Menten microbial uptake and its Arrhenius rate, MAOC saturation and
the CENTURY stabilization factor, arc and mantle degassing, organic P and plant
uptake, the EPIC partition fractions, and the diffusive weathering-front
thickness. The mineral advection-reaction equation is solved in `pde` with an
upwind scheme for both mineral fraction and particle age under a CFL limit.

**PDE solvers (notes 112-128 KB, `src/pde.*`).** The continuum models are
solved numerically in plain C++17, with no PDE library:

- `solve_blatter_pattyn` discretises the first-order membrane/shear balance
  `4 d/dx(mu du/dx) + d/dz(mu du/dz) = rho g ds/dx` on an x-z slice, computes
  Glen viscosity from the strain-rate invariant, applies a stress-free surface,
  no-slip bed and zero-gradient outflow, and iterates with SOR sweeps. For an
  isothermal constant-slope slab it reproduces the SIA profile
  `u = 2A/(n+1) (rho g s')^n H^(n+1)` to within 6% at 61 x 31.
- `solve_stefan_1d` uses an enthalpy formulation with an explicit diffusive
  update and sub-stepping for stability; it tracks the liquid fraction and the
  thaw front, and matches the erfc similarity solution for pure conduction.
- `solve_phase_field_1d` integrates the Allen-Cahn-type damage evolution
  `dphi/dt = -M[-2(1-phi)W+ + Gc/ell(phi - ell^2 lap phi)]` with clamping.
- `evolve_firn_column` integrates the Arthern densification law over a column
  using the depth-integrated overburden pressure.

**Corrected grounding line, KPP and closures.** The mangled grounding-line
expression is replaced by the Schoof (2007) flux
`Q_G = [A (rho_i g)^(n+1)(1-rho_i/rho_w)^n/(4^n ell)]^(ell/(ell+1)) H_G^((ell n+3 ell+1)/(ell+1)) chi^(ell n/(ell+1))`,
with migration `Q_G/H_G` and `dQ_G/dH_G`; at n = ell = 3 this reproduces the
published `H_G^19/4` scaling. KPP now has the cubic shape function
`G_lambda = a1 sigma + a2 sigma^2 + a3 sigma^3` and `V_t = w_lambda G_lambda`;
the `w_lambda` coefficients themselves remain external (Large et al. appendix B).
The empirical laws have explicit sensitivities: `d(pCO2max)/dq`, both stream-power
partials, and `dP_0/dq`. Perirhizal hydraulics now uses the van Genuchten-Mualem
conductivity, its effective saturation and the closed-form matric flux potential,
so the conductance is evaluated without quadrature.

Not reducible to closed form without further data: the full-Stokes equations
(only the first-order Blatter-Pattyn reduction is solved) and SHAKTI's
unspecified creep-closure term.

Two formulas in the notes are dimensionally inconsistent and are implemented
in their standard forms: the mixed-layer depth
`int 2 N^2 dz / int N^2 dz` (dimensionless) is replaced by the `N^2` centroid,
and the internal-wave dispersion's extra `g^2/(rho^2 N^2)` factor is dropped.

**Separating trend from seasonality.** STL (`y = T + S + R`) with cycle-subseries
LOESS smoothing, a low-pass moving-average filter, and bisquare robustness
iterations; monthly climatology and standardised anomalies; Mann-Kendall with
tie correction and the Hamed-Rao autocorrelation variance correction; optional
trend-free prewhitening; Sen's slope with the rank-based confidence interval;
and a within-block seasonal Mann-Kendall test.

## The intent pipeline

`query` implements the architecture from the notes: natural text is parsed into
a structured IR (`process`, `region`, `start_year`, `end_year`), the IR is
validated against the catalog whitelist, and deterministic code selects and
analyses the variables. An external language model can be used purely as the
intent parser: pass its JSON to `query --intent ir.json` or `intent --file ir.json`.
Unknown processes, regions or variable ids are rejected rather than silently
mapped to something else.

## NASA Earthdata ingestion

The C++ core never links a NetCDF, HDF5 or JSON library. Two mechanisms cover real
NASA products:

1. **Classic NetCDF-3 (CDF-1/CDF-2) is parsed by this project.** `src/netcdf.cpp`
implements the header grammar, dimension/variable/attribute records, big-endian
decoding for byte/short/int/float/double, record variables and time units
("days since ...", "months since ...", "years since ..."). Use `data inspect` to
dump a file and `data export` to write the CSV schema directly.
2. **NetCDF-4/HDF5 is converted by a data-side wrapper.** `tools/nasa_fetch.py` uses
`earthaccess` (CMR search and authenticated download), `xarray`/`netCDF4` (HDF5
reads) and `numpy`/`pandas` (area weighting, extent, decimal years), then writes
`time,value` or grid CSV plus a JSON manifest. C++ invokes it with `climate data
search|fetch|convert`; the shell command is quoted and shown with `--print-command`.

```sh
python3 -m venv .venv
./.venv/bin/pip install earthaccess xarray netCDF4
export EARTHDATA_USERNAME=your_username
export EARTHDATA_PASSWORD=your_password

./climate data list
./climate data fetch --dataset sea_ice_extent --start 2024-01-01 --end 2024-12-31 \
        --out data/nasa --max 12
./climate analyze --file data/nasa/sea_ice_extent_monthly.csv --variable sea_ice_extent \
        --units "million km2"
```

`EARTHDATA_USERNAME`/`EARTHDATA_PASSWORD` (or `~/.netrc`) are read only by the
wrapper. `--dry-run` prints the planned query without networking, which is what
runs in CI without credentials.

Dataset keys and their CMR collections live in `src/fetch.hpp`; the short name,
version and variable list can be overridden on the command line, so a renamed
collection never requires a code change. Versions are deliberately left
unpinned because CMR returns nothing when a pinned version no longer exists;
pass `--version` to pin one.

`data/nasa/` contains two small NetCDF fixtures used by the tests and the
examples: `sample_sst4.nc` (NetCDF-4, 180 months of a 45x90 SST field with a
0.02 degC/yr trend) and `sample_sst3.nc` (the same data truncated to 24 months,
written as NETCDF3_CLASSIC), plus `sample_seaice4.nc` for the extent path. They
let the reader and converter be verified without downloading anything.

## Data

`data/` holds synthetic datasets generated by `tools/generate_demo.cpp`:
realistic monthly SST, air temperature, CO2, NDVI and precipitation series with
known trends, gaps, an injected outlier, a sea ice concentration grid, buoy
matchups, NASA Team tie points, precipitation estimates and an AMOC velocity
section. They exist to exercise every code path and to make the tests
verifiable. They are not NASA products.

The default NASA Team tie points and the bundling coefficients are
representative example values. For scientific use, replace them with published
tie points and matchup-derived coefficients from the sensor and period you are
processing: `calibrate` and `retrieve sst --file` regenerate them from data.

## Verification

`make test` exercises the formulas directly: FFT round-trips, LOESS/STL recovery,
harmonic coefficient and derivative recovery, Mann-Kendall/Sen on known series,
NLSST coefficient recovery, NASA Team concentrations, GPCP weighting, spherical
areas, JSON/CSV round-trips, intent parsing, the NetCDF-3 reader and exporter,
and Mie efficiencies cross-checked against an independent implementation to
1e-6. The advanced modules are checked the same way: RK4 convergence to the
two-layer equilibrium, Planck forward/inverse round-trip, two-stream energy
closure, Penman-Monteith closure, closed-form values for the spatial
statistics, CAPE/CIN on a parcel profile, the photolysis integral, the falloff
constant and the oceanic mixing and spectral formulas. `data/nasa/sample_sst3.nc` is parsed and exported by the C++ reader;
`sample_sst4.nc` and `sample_seaice_polar4.nc` are converted by the wrapper
(area-weighted SST, 15% extent on a curvilinear grid). A live CMR search with
the bundled registry returns real NSIDC granules. Downloading requires Earthdata
credentials, which are not stored in this repository.

## Notes and limits

* Time columns are numeric: fractional years (e.g. `2024.9583`) are expected for
  `analyze` and `co2`; other spacings work if `--period` and the window options
  are set consistently.
* `analyze` assumes the requested seasonal period; it does not detect it.
* Missing values are linearly interpolated before decomposition, and the count
  is reported.
* STL requires at least two full seasonal cycles; shorter series fall back to a
  calendar-phase climatology and the report says so.
* The Mann-Kendall p-value uses the normal approximation with continuity
  correction, appropriate for the series lengths this tool targets.
