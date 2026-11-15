#!/usr/bin/env python3
"""Pull real Earth-science series from NASA endpoints.

Nothing here is synthesized: every row written by this tool came from the URL
recorded in the sidecar metadata, and the sidecar keeps the retrieval time so a
cached series can never be mistaken for a fresh one.

Sources
  gistemp      NASA GISS GISTEMP v4 land-ocean monthly temperature anomaly
  power        NASA POWER monthly point API (lat/lon, any of its parameters)
  seaice       NASA NSIDC DAAC Sea Ice Index monthly extent (north or south)
  earthaccess  any Earthdata dataset reachable through the earthaccess library
               (requires Earthdata Login credentials in the environment)

Usage
  python3 tools/nasa_live.py status
  python3 tools/nasa_live.py pull --source power --lat 72 --lon -40 --start 2000 --end 2024
  python3 tools/nasa_live.py pull --source gistemp
  python3 tools/nasa_live.py pull --source seaice --region N
  python3 tools/nasa_live.py pull --source earthaccess --short-name MODISA_L3m_SST

Every pull prints a JSON object describing what was fetched and where it landed.
"""

import argparse
import csv
import datetime as dt
import io
import json
import os
import sys
import urllib.error
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
UA = {"User-Agent": "trend-detective/1.0 (climate workbench; data provenance recorded)"}

SOURCES = {
    "gistemp": {
        "label": "GISTEMP v4 global mean anomaly",
        "provider": "NASA Goddard Institute for Space Studies",
        "url": "https://data.giss.nasa.gov/gistemp/tabledata_v4/GLB.Ts+dSST.csv",
        "variable": "temperature anomaly",
        "units": "degC vs 1951-1980",
        "needs_auth": False,
        "cadence": "monthly, 1880-present",
    },
    "power": {
        "label": "NASA POWER surface meteorology",
        "provider": "NASA Langley Research Center POWER project",
        "url": "https://power.larc.nasa.gov/api/temporal/monthly/point",
        "variable": "T2M",
        "units": "degC",
        "needs_auth": False,
        "cadence": "monthly, 1981-present",
    },
    "seaice": {
        "label": "Sea Ice Index monthly extent",
        "provider": "NASA NSIDC DAAC (G02135)",
        "url": "https://noaadata.apps.nsidc.org/NOAA/G02135/{hemisphere}/monthly/data/",
        "variable": "sea ice extent",
        "units": "million km^2",
        "needs_auth": False,
        "cadence": "monthly, 1978-present",
    },
    "earthaccess": {
        "label": "Earthdata Login granules",
        "provider": "NASA Earthdata (any DAAC)",
        "url": "https://cmr.earthdata.nasa.gov/search",
        "variable": "dataset dependent",
        "units": "dataset dependent",
        "needs_auth": True,
        "cadence": "dataset dependent",
    },
}


def fetch(url, timeout=60):
    request = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", "replace")


def fetch_json(url, timeout=60):
    return json.loads(fetch(url, timeout))


def decimal_year(year, month):
    return round(float(year) + (float(month) - 0.5) / 12.0, 4)


def pull_gistemp():
    text = fetch(SOURCES["gistemp"]["url"])
    rows = []
    reader = csv.reader(io.StringIO(text))
    header_seen = False
    for row in reader:
        if not row:
            continue
        if not header_seen:
            if row[0].strip().lower() == "year":
                header_seen = True
            continue
        try:
            year = int(row[0])
        except (ValueError, IndexError):
            continue
        for month in range(1, 13):
            cell = row[month].strip() if month < len(row) else ""
            cell = cell.replace("***", "").strip()
            if not cell:
                continue
            try:
                rows.append((decimal_year(year, month), float(cell)))
            except ValueError:
                continue
    return rows, {"url": SOURCES["gistemp"]["url"], "variable": "anomaly",
                  "units": SOURCES["gistemp"]["units"]}


def pull_power(lat, lon, start, end, parameter):
    # The POWER monthly endpoint takes calendar years, not ISO dates.
    start_year = start[:4] if start[:4].isdigit() else "1981"
    end_year = end[:4] if end[:4].isdigit() else dt.date.today().strftime("%Y")
    query = ("?parameters=%s&community=RE&longitude=%s&latitude=%s"
             "&start=%s&end=%s&format=CSV" % (parameter, lon, lat, start_year, end_year))
    url = SOURCES["power"]["url"] + query
    text = fetch(url)
    rows = []
    header = None
    months = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
              "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("-END HEADER-"):
            continue
        parts = [p.strip() for p in stripped.split(",")]
        if header is None:
            if "YEAR" in parts:
                header = parts
            continue
        if len(parts) != len(header):
            continue
        record = dict(zip(header, parts))
        if record.get("PARAMETER") != parameter:
            continue
        try:
            year = int(record["YEAR"])
        except (KeyError, ValueError):
            continue
        for offset, name in enumerate(months, start=1):
            raw = record.get(name, "")
            if raw in ("", "-999", "-999.0"):
                continue
            try:
                value = float(raw)
            except ValueError:
                continue
            rows.append((decimal_year(year, offset), value))
    return rows, {"url": url, "variable": parameter,
                  "units": "degC" if parameter == "T2M" else "API native",
                  "latitude": lat, "longitude": lon}


def pull_seaice(region):
    hemisphere = "north" if region.upper() == "N" else "south"
    base = SOURCES["seaice"]["url"].format(hemisphere=hemisphere)
    merged = {}
    for month in range(1, 13):
        name = "%s_%02d_extent_v4.0.csv" % ("N" if hemisphere == "north" else "S", month)
        text = fetch(base + name)
        for line in text.splitlines():
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 5 or not parts[0].isdigit():
                continue
            year, _, _, _, extent = parts[:5]
            try:
                value = float(extent)
            except ValueError:
                continue
            try:
                merged[(int(year), int(parts[1]))] = value
            except ValueError:
                continue
    rows = [(decimal_year(year, month), value)
            for (year, month), value in sorted(merged.items())]
    # Report a URL that was actually fetched instead of the template.
    return rows, {"url": base + ("N_01_extent_v4.0.csv" if hemisphere == "north"
                                 else "S_01_extent_v4.0.csv"),
                  "variable": "extent_%s" % region.upper(),
                  "units": "million km^2"}


def pull_earthaccess(short_name, start, end):
    try:
        import earthaccess  # noqa: F401
        import xarray as xr
    except ImportError as error:
        raise RuntimeError("earthaccess/xarray are not installed in this interpreter: %s"
                           % error)
    try:
        auth = earthaccess.login(strategy="environment", persist=False)
    except Exception as error:
        raise RuntimeError("Earthdata Login is required for this dataset and no "
                           "credentials were found in the environment: %s" % error)
    if auth is not None and getattr(auth, "authenticated", False) is False:
        raise RuntimeError("Earthdata Login rejected the supplied credentials")
    results = earthaccess.search_data(short_name=short_name, temporal=(start, end),
                                      count=8)
    if not results:
        raise RuntimeError("no granules matched %s for %s..%s" % (short_name, start, end))
    files = earthaccess.download(results, local_path=os.path.join(DATA, "earthaccess"))
    if not files:
        raise RuntimeError("search returned granules but none could be downloaded")
    datasets = []
    for path in files:
        try:
            datasets.append(xr.open_dataset(path))
        except Exception:
            continue
    if not datasets:
        raise RuntimeError("downloaded granules could not be opened as NetCDF")
    merged = xr.concat(datasets, dim="time") if len(datasets) > 1 else datasets[0]
    variable = None
    for name, data in merged.data_vars.items():
        if data.ndim >= 2:
            variable = name
            break
    if variable is None:
        variable = list(merged.data_vars)[0]
    series = merged[variable]
    for dim in list(series.dims):
        if dim != "time":
            series = series.mean(dim=dim, skipna=True)
    times = series["time"].values
    rows = []
    for stamp, value in zip(times, series.values):
        try:
            when = dt.datetime.fromisoformat(str(stamp)[:19])
        except ValueError:
            continue
        rows.append((decimal_year(when.year, when.month), float(value)))
    return rows, {"url": "earthaccess search: %s" % short_name,
                  "variable": variable, "units": str(series.attrs.get("units", "")),
                  "files": len(files)}


def whoami():
    """Report whether the Earthdata credentials in the environment can log in.

    A rejected password, a missing token and an uninstalled library are reported
    as different errors so the page can tell the user which one it is.
    """
    try:
        import earthaccess
    except ImportError as error:
        return {"ok": False, "authenticated": False,
                "error": "earthaccess is not installed in this interpreter: %s" % error}
    user = os.environ.get("EARTHDATA_USERNAME")
    token = os.environ.get("EARTHDATA_TOKEN")
    has_netrc = os.path.exists(os.path.expanduser("~/.netrc"))
    if not token and not (user and os.environ.get("EARTHDATA_PASSWORD")) and not has_netrc:
        return {"ok": False, "authenticated": False,
                "error": "no Earthdata credentials were supplied",
                "hint": "Enter your Earthdata username and password, or paste a bearer "
                        "token from https://urs.earthdata.nasa.gov/users/credentials"}
    try:
        auth = earthaccess.login(strategy="environment", persist=False)
    except Exception as error:
        return {"ok": False, "authenticated": False, "error": str(error),
                "hint": "If your Earthdata account uses multi-factor authentication, a "
                        "password cannot be used: paste a bearer token instead."}
    authenticated = bool(auth is not None and getattr(auth, "authenticated", False))
    name = None
    for attribute in ("username", "user_name"):
        value = getattr(auth, attribute, None)
        if value:
            name = str(value)
            break
    return {"ok": authenticated, "authenticated": authenticated,
            "username": name or (user or None),
            "method": "token" if token else ("password" if user else "netrc"),
            "hint": None if authenticated else
                    "Earthdata rejected these credentials. Check the password, or use a "
                    "bearer token if the account has multi-factor authentication."}


def write_series(key, rows, meta):
    rows = [(round(float(t), 4), float(v)) for t, v in rows if v is not None]
    rows.sort()
    if len(rows) < 4:
        raise RuntimeError("only %d usable rows were returned" % len(rows))
    os.makedirs(DATA, exist_ok=True)
    csv_path = os.path.join(DATA, "live_%s.csv" % key)
    meta_path = os.path.join(DATA, "live_%s.json" % key)
    with open(csv_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time", "value"])
        for time_value, value in rows:
            writer.writerow([("%.4f" % time_value).rstrip("0").rstrip("."), value])
    record = dict(meta)
    record.update({
        "key": key,
        "file": "data/live_%s.csv" % key,
        "samples": len(rows),
        "first_time": rows[0][0],
        "last_time": rows[-1][0],
        "first_value": rows[0][1],
        "last_value": rows[-1][1],
        "retrieved_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    })
    with open(meta_path, "w") as handle:
        json.dump(record, handle, indent=2)
    return record


def status():
    cached_series = []
    for name in sorted(os.listdir(DATA)) if os.path.isdir(DATA) else []:
        if not (name.startswith("live_") and name.endswith(".json")):
            continue
        try:
            cached_series.append(json.load(open(os.path.join(DATA, name))))
        except Exception:
            continue
    cached_series.sort(key=lambda item: item.get("retrieved_at", ""))

    out = []
    for key, spec in SOURCES.items():
        prefix = "ea_" if key == "earthaccess" else key
        matches = [item for item in cached_series
                   if str(item.get("key", "")).startswith(prefix)]
        entry = dict(spec)
        entry["key"] = key
        entry["cached"] = matches[-1] if matches else None
        entry["series"] = matches
        out.append(entry)
    credentials = bool(os.environ.get("EARTHDATA_TOKEN") or
                       (os.environ.get("EARTHDATA_USERNAME") and
                        os.environ.get("EARTHDATA_PASSWORD")) or
                       os.path.exists(os.path.expanduser("~/.netrc")))
    return {"sources": out, "earthaccess_credentials": credentials,
            "cached_series": cached_series,
            "interpreter": sys.executable, "python": sys.version.split()[0]}


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("status")
    sub.add_parser("whoami")
    pull = sub.add_parser("pull")
    pull.add_argument("--source", required=True, choices=sorted(SOURCES))
    pull.add_argument("--lat", type=float, default=72.0)
    pull.add_argument("--lon", type=float, default=-40.0)
    pull.add_argument("--start", default="2000-01-01")
    pull.add_argument("--end", default=dt.date.today().isoformat())
    pull.add_argument("--parameter", default="T2M")
    pull.add_argument("--region", default="N")
    pull.add_argument("--short-name", default="")
    args = parser.parse_args(argv)

    if args.command == "whoami":
        print(json.dumps(whoami(), indent=2))
        return 0

    if args.command == "status":
        print(json.dumps(status(), indent=2))
        return 0

    try:
        if args.source == "gistemp":
            rows, meta = pull_gistemp()
            key = "gistemp"
        elif args.source == "power":
            rows, meta = pull_power(args.lat, args.lon, args.start, args.end, args.parameter)
            key = "power_%s_%s_%s" % (args.parameter.lower(),
                                      ("%g" % args.lat).replace("-", "m"),
                                      ("%g" % args.lon).replace("-", "m"))
        elif args.source == "seaice":
            rows, meta = pull_seaice(args.region)
            key = "seaice_%s" % args.region.lower()
        else:
            if not args.short_name:
                raise RuntimeError("--short-name is required for the earthaccess source")
            rows, meta = pull_earthaccess(args.short_name, args.start, args.end)
            key = "ea_%s" % args.short_name.lower()
        meta["source"] = args.source
        record = write_series(key, rows, meta)
    except (urllib.error.URLError, RuntimeError, OSError, ValueError) as error:
        print(json.dumps({"ok": False, "source": args.source, "error": str(error)}))
        return 2
    print(json.dumps({"ok": True, "source": args.source, "series": record}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
