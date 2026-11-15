# Complete guide

Everything in the software, screen by screen and endpoint by endpoint: what it is, what
each control and number means, and how it serves the challenge — *select variables that
indicate change, visualize them over time, and determine what is changing, where, how
much, and whether it is significant.*

---

## 0. Running it

```bash
make -j4                              # build the engine and the server
./climate ui --port 3002 --host 0.0.0.0   # serve the workbench
```

The terminal prints the address. Open it; there is no login, no build step and no
database. Live pulls need outbound HTTPS. Earthdata datasets additionally need a free
NASA Earthdata Login (the Fetch tab can verify it before any download).

Other entry points that share the same engine:

```bash
./climate analyze --file data/sst_monthly.csv     # trend report in the terminal
./climate --help                                  # every command
make test                                         # self-checking engine tests
make formulas                                     # regenerate the formula registry
```

---

## 1. The five challenge tabs

### 1.1 Fetch data — *connect, search, pull*

**Purpose:** get a real measured series into the workbench.

- **Pull parameters** — latitude, longitude, POWER parameter, start and end dates, sea-ice
  hemisphere, Earthdata short name. These apply to the source cards below.
- **Source cards** — one per source. **Pull from …** performs a real HTTP request through
  the C++ server and the Python adapter, then reports: the source URL, the sample count,
  the first and last time stamp, units, the file written and the elapsed time, with the
  series charted underneath. A live timer runs while the request is in flight and the
  request aborts at 180 s with an explicit error panel instead of hanging.
  - `gistemp` — NASA GISS global land/ocean temperature anomaly, monthly, 1880 → present.
  - `power` — NASA POWER monthly point API. Any of T2M, T2M_MAX, T2M_MIN, TS,
    PRECTOTCORR, WS2M, RH2M, ALLSKY_SFC_SW_DWN at any latitude/longitude. No login.
  - `seaice` — NSIDC Sea Ice Index monthly extent, north or south, 1978 → present.
  - `earthaccess` — any Earthdata dataset by CMR short name (default `MODISA_L3m_SST`),
    opened with `earthaccess` + `xarray` + `netCDF4`, averaged over space to a monthly
    series arithmetic mean.
- **Earthdata Login** — optional, per request. **Verify Earthdata login** calls the
  `whoami` path and distinguishes three outcomes: no credentials supplied, credentials
  rejected by NASA, and an account that cannot use a password (multi-factor accounts must
  paste a bearer token from `urs.earthdata.nasa.gov`). Credentials are kept in the tab,
  sent with the single pull, handed to the child process through its environment, and
  never written to disk, logs or the page.
- **NASA product planner** — choose a product, a date range and a granule limit; **Generate
  fetch plan** renders the exact command and a dry-run manifest inline (exit code
  included), and **Query CMR API** performs a live metadata search that downloads nothing.

**Feeds the challenge:** this is the "variables measured by Earth science missions or
produced by Earth science models" part, with provenance attached.

### 1.2 Select variables — *what kind of change does this measure?*

Sixteen variables: sea-surface temperature, surface air temperature, CO₂, CH₄, NDVI,
sea-ice extent, sea-ice concentration, precipitation, ocean heat content, TOA energy
imbalance, snow cover, ozone, aerosol optical depth, soil moisture, ocean colour, sea
level.

Each card shows the identifier, the name and units, **the kind of change it measures**
(e.g. *"Trend in sea-ice extent (shrinking / expanding)"*), its processes and its regions.
**Analyse this variable →** opens the analysis tab, loads that variable's record, confirms
the data gate and runs the whole chain — acquisition, seasonality removal, trend and
significance, report — filling the four-question bar. **Plot over time** jumps to the
series view.

**Feeds the challenge:** this is the "difficult to select variables that could indicate
change" part, answered by labelling the change each variable carries.

### 1.3 Visualize over time

One card per series the workbench can reach (bundled fixtures and live pulls): the real
values as a line, the mean, the range, the Sen slope and the end-to-end percentage change.
Below them, a **spatial trend field** — a 12 × 24 demonstrated grid showing per-cell Sen
slope, colour-coded hot/cold, with Getis-Ord Gi* hotspots counted. It is labelled as a
synthetic demonstration of the grid mechanics, not a downloaded product.

**Feeds the challenge:** the visual answer to "how it changes over time", and the only
place where a spatial field is drawn.

### 1.4 What / where / how much / significant

The centre of the workbench.

- **Four-question bar** — filled by every analysis:
  - **What is changing** — the variable and the direction Mann–Kendall found.
  - **Where is it changing** — the region, or the grid support of the record.
  - **How much** — Sen's slope per step and per decade, in the variable's units.
  - **Is it significant** — yes/no with the p-value.
- **Regional contrast · one variable, two places** — choose a POWER parameter, two
  coordinate pairs and a year range. Two real POWER series are pulled **in parallel**,
  each gets its own Sen slope, seasonal Mann–Kendall test, confidence interval, total and
  relative change, and a sparkline. The verdict states whether the two regions are
  opposed, whether both are significant, or whether one cannot be separated from no
  trend — the exact situation the challenge describes. On failure it explains what
  blocked it and offers **Retry**, **Open Fetch data** and **Open Significance**.
- **Guided workflow** — five steps on one record: 01 data acquisition with coverage,
  gaps and missing values; 02 seasonality removal (STL, or a climatology fallback when
  the record is shorter than two cycles); 03 trend and significance with the estimator
  ensemble (Sen, OLS, Kendall τ, and others, scored by residual RMSE, significance and
  sign agreement); 04 spatial analysis; 05 final report with WHAT / WHERE / HOW MUCH /
  SIGNIFICANT in words, plus print, copy-JSON and CSV export.
- **Regional-pattern checks** — every estimator sees the same data; disagreement is shown
  rather than hidden.

### 1.5 Significance

Seven tests on a chosen variable, each labelled with the question it answers:

| Question | Test | What a "yes" means |
|---|---|---|
| Is there a monotonic trend at all? | Mann–Kendall on the deseasonalised series | S, z and p say the ranks climb or fall consistently |
| Does the raw record already show it? | Mann–Kendall without removing the trend first | The signal is not an artefact of preprocessing |
| Is the seasonal cycle hiding it? | Seasonal Mann–Kendall, month by month | The change survives comparing like months |
| Is autocorrelation inflating the test? | Hamed–Rao variance correction | The variance inflation factor shows how much the test was penalised |
| How large is the change, with uncertainty? | Sen's slope and its rank-based 95 % interval | The interval excludes zero |
| Do the estimators agree on the sign? | Accuracy-weighted ensemble | Every estimator points the same way |
| Could a trend this size be detected at all? | Minimum detectable trend, 80 % power | The record is long and clean enough to resolve a change of that size |

**Feeds the challenge:** the "significant" part, made falsifiable — and the last row
prevents the classic mistake of reading "not significant" as "no change".

---

## 2. Deep dives

- **Time series lab** — paste or upload a CSV, choose time and value columns, and run
  decomposition and trend tests directly. Includes an upload control with an English file
  label, so the choice cannot be blocked by browser locale.
- **CO₂ curve fit** — the Thoning-style polynomial-plus-harmonic fit used for the Keeling
  curve: long-term trend, seasonal cycle, growth rate from the analytic derivative, with
  residual diagnostics.
- **Retrievals** — single-footprint physics: NLSST sea-surface temperature from brightness
  temperatures with a first guess, zenith angle and mirror-side terms; NDVI / EVI / SAVI
  from reflectances; NASA Team sea-ice concentration from passive-microwave channels;
  GPI precipitation from cold-cloud duration. Below them, **How the retrieved variable
  changes over time** plots the full record of the same quantity (SST, NDVI, sea ice,
  precipitation, or a live GISTEMP/Arctic-ice pull) with STL trend, Sen slope and a
  Mann–Kendall verdict — so a retrieval is never a dead-end number.
- **Formula calculator** — the physics core: energy balance, radiative forcing, Mie
  scattering, wind-driven transport, geostrophy, vegetation indices, selected in any
  combination with per-formula inputs.
- **Information theory** — entropy, mutual information, transfer entropy between series,
  for asking whether two variables carry the same signal.
- **Uncertainty** — propagation through the retrieval chain, error bars on derived
  quantities, and the effect of uncertainty on a trend verdict.

---

## 3. Engine tools

- **Command console** — every CLI command and diagnostic executed by the same binary, one
  at a time or as a batch, with exit codes.
- **Calculation modes** — the 944 public functions, grouped by the question they answer:
  01 what is changing (observation, retrieval, Earth-system variables), 02 how much
  (trend estimators, change detection, cycles), 03 where (space, region, coherence,
  attribution), 04 is it significant (uncertainty, power, validation, ensembles), 05 why
  (physics, dynamics, transport), 06 numerics. Category chips filter to one group, search
  filters by name, and each group expands to its function signatures. File names are not
  shown: the grouping is by question, not by source file.
- **All calculations** — the complete lab. 944 functions are bound with a generated call
  site, argument conversion and result formatter; pick any subset, give inputs, run them
  as one batch. Default inputs execute 749 of them; the remaining 195 reject their
  defaults by design (domain guards such as negative radii), and the rejection is
  reported per formula instead of crashing the batch.
- **System** — module overview, catalog counts, build facts and the endpoint list.

---

## 4. Regional contrast: why it exists and how to read it

The challenge states the core difficulty: one process can drive opposite signs in two
regions. A global average cannot show this; a single time series cannot show this. The
contrast panel pulls the *same* variable at two coordinates and runs an independent
pipeline for each:

1. monthly values from POWER (no login),
2. Sen's slope in units per year with its 95 % interval,
3. seasonal Mann–Kendall p-value (months compared against the same month in other years),
4. total change over the record and its relative size,
5. a verdict that names the situation: opposite and both significant, opposite with one
   insignificant, same direction, one significant, or neither.

Worked example from real data (T2M, 1990–2024, 72 °N / 40 °W versus 70 °S / 0 °E):
Arctic +0.041 °C/yr with p = 0.043; Antarctic −0.072 °C/yr with p < 0.0001. A global mean
of this pair would report a mild warming and delete both stories.

---

## 5. Endpoints

| Endpoint | Method | Purpose |
|---|---|---|
| `/` | GET | The workbench page |
| `/api/analyze` | POST | Full series analysis: decomposition, MK, Sen, ensemble, report |
| `/api/co2` | POST | CO₂ curve fit with growth rate |
| `/api/retrieve/sst` | POST | NLSST retrieval from brightness temperatures |
| `/api/retrieve/index` | POST | NDVI / EVI / SAVI from reflectances |
| `/api/retrieve/seaice` | POST | NASA Team sea-ice concentration |
| `/api/retrieve/precip` | POST | GPI precipitation |
| `/api/physics` | POST | Physics formula evaluation |
| `/api/lab` | GET | The 944-function registry, grouped |
| `/api/lab/run` | POST | Batch execution with per-formula results |
| `/api/catalog` | GET | Variables with the kind of change each measures |
| `/api/visualize` | GET | Every series with values, slope and percent change |
| `/api/grid-trend` | GET | Per-cell Sen slope and Gi* hotspots (synthetic demo field) |
| `/api/contrast` | POST | Two-region live POWER comparison and sign partition |
| `/api/live/status` | GET | Sources, cached pulls and credential state |
| `/api/live/pull` | POST | Real pull (gistemp, power, seaice, earthaccess) |
| `/api/earthdata/check` | POST | Verify Earthdata login without downloading |
| `/api/nasa/datasets` | GET | NASA product registry |
| `/api/nasa/plan` | POST | Dry-run fetch plan and exact command |
| `/api/nasa/search` | POST | Live CMR metadata search |
| `/api/data/file` | GET | One cached CSV by name |
| `/api/overview`, `/api/help`, `/api/formulas`, `/api/query`, `/api/run`, `/api/info`, `/api/sstunc` | varied | System overview, CLI help, registry, intent query, CLI runner, information theory, SST uncertainty |

---

## 6. Data on disk

- `data/live_<key>.csv` and `.json` — every live pull, with its source URL, sample count,
  time range and retrieval timestamp. The JSON is the provenance record.
- `data/*_monthly.csv` — bundled fixture series used by the demo analysis and tests.
- `data/nasa/` — dry-run manifests and any fetched NASA product output.
- `data/nasa_team_model.csv`, `nlsst_model.csv` — retrieval coefficients regressed from
  the bundled matchup files.

Nothing in `data/` is regenerated silently: a pull overwrites only its own key.

---

## 7. Tests

```bash
make test                                             # self-checking engine suite
.venv-test/bin/cmake -S tests -B .build/tests && .venv-test/bin/cmake --build .build/tests -j4
cd .build/tests && ../../.venv-test/bin/ctest --output-on-failure
.venv-test/bin/pytest tests/properties tests/bdd      # Hypothesis + BDD
```

- `tests/test_all.cpp` — self-checking assertions across the engine, including the formula
  lab report (944 bound, 749 executed, 195 rejected defaults).
- `tests/gtest_suite.cpp` — 17 GoogleTest cases over trends, retrievals and I/O.
- `tests/rapidcheck_suite.cpp` — properties such as "a constant series has no trend" and
  "EVI stays inside the physical bounds".
- `tests/fp_diff.cpp` — floating-point differential and determinism checks.
- `tests/fuzz_lab.cpp` — libFuzzer entry point for the lab runner (needs clang).
- `tests/properties/test_workbench_api.py` — Hypothesis property tests against the API.
- `tests/features/` + `tests/bdd/` — Gherkin scenarios for "change is always occurring",
  executed with pytest-bdd.

---

## 8. Reading the numbers

- **Sen slope** — the median of all pairwise slopes. Robust to outliers and to a handful
  of bad months; reported per step (a month in monthly data) and per decade.
- **Mann–Kendall S, τ, z, p** — a rank-based monotonic-trend test. No assumption of
  normality; ties are corrected; Hamed–Rao inflates the variance when the series is
  autocorrelated.
- **Seasonal Mann–Kendall** — compares each month against the same month in other years,
  so the seasonal cycle cannot masquerade as a trend.
- **95 % interval on the Sen slope** — if it excludes zero, the direction is resolved.
- **Estimator ensemble** — several estimators on the same data, scored; agreement is
  printed so a single method cannot hide disagreement.
- **Minimum detectable trend** — 2.802 · σ · √(12/n) / span for 80 % power at α = 0.05.
  A verdict of "no significant trend" alongside a large minimum detectable trend means the
  record cannot answer the question, not that nothing happened.
- **Relative change** — total change divided by |mean|, so the size is readable across
  variables with different units.

---

## 9. Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| "Calling…" then a timeout panel | endpoint unreachable from this machine | check outbound HTTPS, retry, or use an open source (GISTEMP/POWER/NSIDC need no login) |
| Earthdata "no credentials were supplied" | empty form and no exported variables | type credentials, or paste a bearer token, then **Verify Earthdata login** |
| Earthdata "Invalid user credentials" | wrong password, or an MFA account | reset the password, or use a token |
| Contrast panel reports a failure | POWER unreachable, or a latitude with no data | retry, or move the point; error text names the region that failed |
| "No bundled demo for this variable" | that variable has no local fixture | pull it from Fetch data, or upload a CSV in the time-series lab |
| Page loads but no data | server not running | start `./climate ui` in a terminal and keep it open |

---

## 10. Honest limits

- The 12 × 24 spatial trend field is a synthetic demonstration of per-cell Sen slope and
  Gi* hotspots, labelled as such. Every other chart plots pulled or bundled values.
- Earthdata downloads need the user's own credentials and are bounded by NASA's latency.
- POWER's monthly endpoint returns calendar years, so a POWER pull snaps the requested
  dates to year boundaries; that is stated in the returned metadata.
- The database-free design means the workbench is single-user: pulls overwrite their own
  cache files rather than being versioned.
