# Trend Detective — project file

**One sentence:** a workbench that pulls real Earth-science variables, shows how they
change over time, and answers — with numbers and a significance verdict — *what* is
changing, *where*, *how much*, and whether the change is statistically real.

---

## 1. What it does and how it works

The challenge states that change is always occurring in the Earth's interconnected
environmental system, that a measured variable can trend one way in one region and the
opposite way in another even when one process drives both, and that the hard part is
choosing variables that indicate change and then analysing them. The software is built
as the five steps of that sentence.

| Step | Screen | What happens |
|---|---|---|
| 1 | **Fetch data** | Real HTTP pulls from NASA GISTEMP (global land/ocean anomaly), the NASA POWER monthly point API (any parameter, any latitude/longitude), the NSIDC Sea Ice Index, and any Earthdata dataset through `earthaccess`. Every pull writes a CSV plus a JSON record with the source URL, sample count, time range and retrieval time. A NASA product planner returns the exact fetch command and a dry-run manifest before anything is downloaded. |
| 2 | **Select variables** | Sixteen mission/model variables, each card stating what kind of change it measures ("Trend in mean ocean skin temperature", "Trend in sea-ice extent (shrinking / expanding)", "Trend in soil moisture (wetting / drying)"…), its processes, its regions and its units. Clicking a card carries that variable into the analysis and runs it. |
| 3 | **Visualize over time** | Every series the workbench can reach, drawn from its real values: raw record, STL trend, Sen slope, end-to-end percentage change, and a 12×24 spatial trend field for the grid path. |
| 4 | **What / where / how much / significant** | The four answers on one screen. A four-cell bar reports the variable and direction, the region, the slope per step and per decade, and the Mann–Kendall verdict. Below it, **Regional contrast** pulls the same variable live from POWER at two coordinates and reports two Sen slopes and two seasonal Mann–Kendall tests side by side — so an opposition between regions is measured instead of averaged away. The guided workflow underneath walks the same record through acquisition, seasonality removal, trend + significance, spatial analysis and a final report. |
| 5 | **Significance** | Seven labelled tests on the chosen record: Mann–Kendall on the deseasonalised series, Mann–Kendall on the raw record, seasonal Mann–Kendall, the Hamed–Rao autocorrelation correction, the Sen slope with its rank-based 95 % interval, the accuracy-weighted estimator ensemble, and the minimum detectable trend for a record of this length. Each row states the question it answers. |

**Engine.** A single C++17 binary (`./climate`) serves the page and does the science:
POSIX sockets and threads for HTTP, an in-tree JSON parser, CSV reader, STL
decomposition, Mann–Kendall with tie correction and Hamed–Rao variance correction,
seasonal Mann–Kendall, Sen slope with rank confidence interval, an accuracy-weighted
trend-estimator ensemble, Theil–Sen-based prewhitening, Getis-Ord Gi* hotspots,
Benjamini–Hochberg false-discovery control, change-point detection, information-theory
measures, uncertainty propagation and a formula registry of 944 public functions.
Python is used only as a data adapter: HTTP fetching, `earthaccess`/`xarray`/`netCDF4`
for Earthdata granules, and nothing else.

**Data path.** Browser → JSON over HTTP → C++ endpoint. The C++ side either computes
directly (analysis, trends, significance, formulas) or runs the Python adapter as a
child process for network access, then parses the returned JSON and the cached CSV.

---

## 2. Benefits

- **A judge or scientist can start from data, not from a formula list.** The first tab
  produces a real series in about a second; the last tab defends it.
- **The four questions are answered explicitly, not implied.** Every analysis ends with
  four filled cells, and the report repeats them in words.
- **Regional opposition is first-class.** Two live POWER series, two independent
  verdicts, one comparison. A single global trend can never hide a falling region.
- **Significance is not a p-value in isolation.** Autocorrelation, the seasonal cycle,
  estimator disagreement and record length are all checked, and the minimum detectable
  trend says whether "no trend" means "no change" or "not enough record".
- **Everything is auditable.** Each pulled series keeps its source URL and retrieval
  time; each formula is callable with your own inputs; each test shows its statistic.
- **Runs offline except for the network calls.** No cloud service, no account, no
  database, no build step for the user.

## 3. Intended impact

For the challenge: a working instrument that turns "select variables that could indicate
change, then analyse them" into five clicks, and that demonstrates the stated physical
fact — same process, opposite regional signs — with real measurements rather than an
assertion.

Beyond the challenge: the same loop (pull → label the kind of change → plot → partition
by region → test) is how operational climate monitoring is actually done, and the tool
keeps every step visible instead of collapsing it into one number.

## 4. Tools, languages, hardware and software

**Languages:** C++17 (all science and serving), JavaScript (the page, no framework),
Python 3 (data adapters and property tests), CMake/CTest, shell, a small amount of SQL-free
CSV/JSON.

**Libraries:** C++ standard library only for the engine — no Eigen, no Boost, no
framework. Python: `urllib`, `csv`, `json` for open endpoints; `earthaccess`, `xarray`,
`netCDF4` for Earthdata. Tests: GoogleTest, RapidCheck, Hypothesis, pytest, pytest-bdd,
gcovr, and a hand-written floating-point differential harness plus a libFuzzer target.

**Tools:** g++, make, CMake/CTest, git, a POSIX shell, a browser for the UI.

**Hardware:** any x86-64 or ARM machine with a POSIX shell. The build is a single
`make -j4`; the server needs no GPU and no more than a few hundred MB of RAM. The demo
was developed and run on a Linux cloud VM.

**Software requirements:** a C++17 compiler, `make`, and Python 3 for live pulls. The
Earthdata path additionally needs the three Python packages above and a free NASA
Earthdata Login (username/password or a bearer token).

## 5. How the project is creative

- **The UI is organised by the question, not by the method.** Tabs are
  Fetch → Select → Visualize → What/where/how much → Significance; the formula list is
  grouped under the same five questions.
- **The palette is a physical argument.** Warm orange carries the accent and the data,
  beige carries prose, brown is the shell; the look is a CRT terminal rendered at modern
  resolution — scanlines, vignette, phosphor glow, hexagon step badges — instead of a
  generic dashboard.
- **The contrast panel is a demonstration, not a chart.** It exists to prove the
  challenge's own premise on demand, in the region pair the user chooses.
- **Honesty is designed in.** Synthetic fixtures are labelled as demonstrations, live
  pulls carry their URL and time, and "no significant trend" is reported together with
  the smallest trend the record could have detected.
- **One binary does the science.** 944 public functions are generated into a callable
  registry, so the page, the CLI and the tests exercise exactly the same code.

## 6. Factors considered

- **Scientific validity:** non-parametric tests (Mann–Kendall, Sen) instead of least
  squares; tie correction; Hamed–Rao variance correction for serial correlation;
  seasonal testing for monthly data; confidence intervals printed with every slope;
  minimum detectable trend printed with every verdict.
- **The stated regional-opposition problem:** solved by design, with two regions
  measured independently.
- **Data provenance:** every series records where it came from and when.
- **Credentials:** Earthdata credentials are used per request, passed to the child
  process through its environment, never written to disk, logs or the page, and can be
  verified before a download.
- **Failure modes:** network calls have visible timers, 180-second aborts and explicit
  error panels; a dead endpoint cannot leave the page saying "Calling…".
- **Accessibility and legibility:** WCAG-audited palette (text ≥ 7:1, graphics ≥ 3:1),
  square geometry, 9–13 px type, keyboard-focusable controls.
- **Reproducibility:** deterministic demo fixtures, fixed seeds in property tests,
  generated formula registry, single-command build and test.
- **Limits, stated plainly:** the spatial grid field is a labelled synthetic
  demonstration of the grid mechanics; Earthdata downloads require the user's own
  credentials; long POWER or Earthdata requests are limited by NASA's own latency.
