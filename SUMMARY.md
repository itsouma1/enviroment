# summary

**Trend Detective** is a climate-change workbench. It pulls real Earth-science
variables, plots how they change over time, and answers four questions with numbers:
what is changing, where it is changing, how much it is changing, and whether the change
is statistically significant.

**The problem it solves.** The challenge is not a lack of data — it is that Earth
science produces too many variables, that any single variable can trend one way in one
region and the opposite way in another, and that "significant" is a statistical claim,
not a visual impression. A tool that only plots a global average actively hides the
regional opposition the challenge describes.

**How it works.** A single C++17 binary serves a five-step page:

1. **Fetch data** — live HTTP pulls from NASA GISTEMP, the NASA POWER monthly API, the
   NSIDC Sea Ice Index, and any Earthdata dataset through `earthaccess`. Every pull
   writes a CSV and a JSON record with the source URL, sample count, time range and
   retrieval time. A NASA product planner prints the exact fetch command and a dry-run
   manifest before anything downloads.
2. **Select variables** — sixteen mission/model variables, each card stating what kind
   of change it measures. One click carries the variable into the analysis.
3. **Visualize over time** — the raw record, the STL trend, Sen's slope and the change
   over the record, drawn from real values.
4. **What / where / how much / significant** — four filled cells plus a **Regional
   contrast** panel that pulls the same variable from POWER at two coordinates and
   reports two Sen slopes and two seasonal Mann–Kendall verdicts side by side.
5. **Significance** — seven labelled tests: Mann–Kendall (deseasonalised and raw),
   seasonal Mann–Kendall, the Hamed–Rao autocorrelation correction, Sen's slope with its
   rank-based 95 % interval, the estimator ensemble, and the minimum detectable trend
   for that record length.

**Evidence it works.** Verified live in the browser during development: GISTEMP 1762
monthly samples (1880 → 2026) in 0.8 s; POWER T2M at 72 °N, 40 °W, 300 samples in 0.7 s;
NSIDC Arctic sea ice, 576 samples (1978 → 2026). The regional contrast at 72 °N versus
70 °S returns +0.041 °C/yr (p = 0.043) against −0.072 °C/yr (p < 0.0001) — an
opposite-signed, both-significant pair on real data, reported as a sign partition rather
than a global mean.

**Stack.** C++17 standard library only for the engine and server (POSIX sockets and
threads, in-tree JSON and CSV), JavaScript for the page with no framework, Python only as
a network/`earthaccess` adapter. 944 public functions are bound into a callable registry
so the page, the CLI and the tests run the same code. Tests: a self-checking C++ binary,
GoogleTest, RapidCheck properties, a floating-point differential harness, a libFuzzer
target, Hypothesis property tests, pytest-bdd scenarios and CMake/CTest orchestration
with gcovr coverage.

**Run it.**

```
make -j4
./climate ui --port 3002 --host 0.0.0.0
```

Open the printed address; the page needs no build step and no account. Earthdata
datasets need a free NASA Earthdata Login, which can be verified from the Fetch tab
(password, or a bearer token for accounts with multi-factor authentication).

**Limits, stated honestly.** The 12×24 spatial trend field is a labelled synthetic
demonstration of the grid mechanics; everything else plots pulled or bundled real values.
Earthdata downloads require the user's own credentials and are bounded by NASA's latency.
