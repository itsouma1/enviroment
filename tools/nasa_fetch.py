#!/usr/bin/env python3
"""NASA Earthdata acquisition and conversion wrapper for the C++ climate tool.

This script is the only component that uses third-party libraries:
  * earthaccess - CMR search and authenticated download
  * xarray / netCDF4 - reading NetCDF-4 and HDF5 granules
  * numpy / pandas - spatial aggregation and time decoding

Everything it produces is plain CSV in the schema the C++ analyser reads, so the
C++ side never depends on any of these libraries.

Credentials are read from the environment (EARTHDATA_USERNAME /
EARTHDATA_PASSWORD) or from ~/.netrc for urs.earthdata.nasa.gov.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

WRAPPER_VERSION = "1.0"


def fail(message: str, code: int = 2) -> None:
    print(f"error: {message}", file=sys.stderr)
    sys.exit(code)


def lazy_import(module: str, hint: str):
    try:
        return __import__(module)
    except ImportError:
        fail(f"{module} is not installed. {hint}")


# ---------------------------------------------------------------------------
# CMR search and download
# ---------------------------------------------------------------------------

def build_search_kwargs(args) -> dict:
    kwargs = {"short_name": args.short_name}
    if args.version:
        kwargs["version"] = args.version
    if args.concept_id:
        kwargs["concept_id"] = args.concept_id
    if args.start and args.end:
        kwargs["temporal"] = (args.start, args.end)
    elif args.start:
        kwargs["temporal"] = (args.start, None)
    if args.bbox:
        try:
            box = tuple(float(value) for value in args.bbox.split(","))
        except ValueError:
            fail("--bbox must be four comma-separated numbers: west,south,east,north")
        if len(box) != 4:
            fail("--bbox must have four values: west,south,east,north")
        kwargs["bounding_box"] = box
    return kwargs


def perform_search(args) -> list:
    earthaccess = lazy_import(
        "earthaccess",
        "Install it with: python3 -m pip install earthaccess (or use the bundled .venv).",
    )
    kwargs = build_search_kwargs(args)
    try:
        granules = earthaccess.search_data(**kwargs)
    except Exception as error:  # network or CMR error
        fail(f"CMR search failed: {error}")
    if args.max:
        granules = granules[: args.max]
    return granules


def granule_summary(granule) -> dict:
    umm = granule.get("umm", {}) if hasattr(granule, "get") else {}
    meta = granule.get("meta", {}) if hasattr(granule, "get") else {}
    size = getattr(granule, "size", None)
    try:
        links = len(granule.data_links())
    except Exception:
        links = None
    return {
        "concept_id": meta.get("concept-id", ""),
        "title": umm.get("GranuleUR", ""),
        "size_mb": round(size / 1e6, 3) if isinstance(size, (int, float)) else None,
        "data_links": links,
    }


def login():
    earthaccess = lazy_import("earthaccess", "Install earthaccess first.")
    strategies = []
    if os.environ.get("EARTHDATA_USERNAME") and os.environ.get("EARTHDATA_PASSWORD"):
        strategies.append("environment")
    strategies.append("netrc")
    last_error = None
    for strategy in strategies:
        try:
            authenticated = earthaccess.login(strategy=strategy)
        except Exception as error:
            last_error = error
            continue
        if authenticated:
            return authenticated
    fail(
        "NASA Earthdata login failed"
        + (f" ({last_error})" if last_error else "")
        + ". Set EARTHDATA_USERNAME and EARTHDATA_PASSWORD, or add an entry for "
        "urs.earthdata.nasa.gov to ~/.netrc."
    )


# ---------------------------------------------------------------------------
# NetCDF -> CSV conversion
# ---------------------------------------------------------------------------

def pick_variable(dataset, requested: str) -> str:
    candidates = [name.strip() for name in (requested or "").split(",") if name.strip()]
    for name in candidates:
        if name in dataset.variables:
            return name
    # Fall back to the largest numeric variable that is not a coordinate.
    best_name = None
    best_score = -1
    for name, variable in dataset.variables.items():
        if name in dataset.coords:
            continue
        if variable.dtype.kind not in "fiu":
            continue
        score = variable.ndim * 1000 + variable.size
        if score > best_score:
            best_score = score
            best_name = name
    if best_name is None:
        fail("no numeric variable found; pass --variable with an explicit name")
    return best_name


def find_dimension(dataset, names) -> str | None:
    for name in names:
        if name in dataset.dims:
            return name
    return None


def cell_areas(latitude, longitude, lat_name, lon_name):
    """Spherical cell areas in km^2 for a regular 1-D latitude/longitude grid."""
    numpy = lazy_import("numpy", "numpy is required for spatial aggregation.")
    import xarray as xr

    lat = numpy.asarray(latitude, dtype=float)
    lon = numpy.asarray(longitude, dtype=float)
    if lat.size < 2 or lon.size < 2:
        fail("grid aggregation needs at least two latitude and two longitude cells")

    def edges(values):
        values = numpy.asarray(values, dtype=float)
        mid = 0.5 * (values[:-1] + values[1:])
        return numpy.concatenate(
            ([values[0] - (mid[0] - values[0])], mid,
             [values[-1] + (values[-1] - mid[-1])])
        )

    lat_edges = edges(lat)
    lon_edges = edges(lon)
    dlat = numpy.abs(numpy.diff(numpy.sin(numpy.deg2rad(lat_edges))))
    dlon = numpy.abs(numpy.deg2rad(numpy.diff(lon_edges)))
    radius = 6371.0088
    area = radius * radius * dlat[:, None] * dlon[None, :]
    return xr.DataArray(area, dims=(lat_name, lon_name))


def cell_areas_2d(latitude, longitude, dims):
    """Spherical cell areas for curvilinear (2-D) latitude/longitude arrays."""
    numpy = lazy_import("numpy", "numpy is required for spatial aggregation.")
    import xarray as xr

    lat = numpy.deg2rad(numpy.asarray(latitude, dtype=float))
    lon = numpy.deg2rad(numpy.asarray(longitude, dtype=float))
    dlat_dy, dlat_dx = numpy.gradient(lat)
    dlon_dy, dlon_dx = numpy.gradient(lon)
    jacobian = numpy.abs(dlat_dy * dlon_dx - dlat_dx * dlon_dy)
    radius = 6371.0088
    area = radius * radius * jacobian * numpy.cos(lat)
    return xr.DataArray(area, dims=dims)


def decimal_years(values):
    numpy = lazy_import("numpy", "numpy is required for time decoding.")
    import pandas as pd

    index = pd.DatetimeIndex(numpy.asarray(values))
    if len(index) == 0:
        return numpy.array([], dtype=float)
    years = index.year.values
    day_of_year = index.day_of_year.values
    leap = ((years % 4 == 0) & (years % 100 != 0)) | (years % 400 == 0)
    days_in_year = numpy.where(leap, 366.0, 365.0)
    return years + (day_of_year - 0.5) / days_in_year


def load_dataset(files):
    xarray = lazy_import("xarray", "Install xarray and netCDF4, or use the bundled .venv.")
    if len(files) == 1:
        return xarray.open_dataset(files[0], decode_times=True, mask_and_scale=True)
    try:
        return xarray.open_mfdataset(
            files, combine="by_coords", decode_times=True, mask_and_scale=True
        )
    except Exception:
        # Granules with duplicate or non-monotonic coordinates (for example OCO-2
        # soundings) are concatenated on a synthetic record dimension instead.
        datasets = [xarray.open_dataset(f, decode_times=True, mask_and_scale=True) for f in files]
        return xarray.concat(datasets, dim="record", combine_attrs="override")


def convert_files(files, out_csv, variable, mode, units, threshold, extra_notes=None):
    numpy = lazy_import("numpy", "numpy is required for conversion.")
    xarray = lazy_import("xarray", "Install xarray and netCDF4, or use the bundled .venv.")

    if not files:
        fail("no input files to convert")

    dataset = load_dataset(files)
    name = pick_variable(dataset, variable)
    field = dataset[name]

    time_name = None
    for candidate in ("time", "Time", "TIME", "date"):
        if candidate in field.dims:
            time_name = candidate
            break
    if time_name is None:
        time_name = next((d for d in field.dims if "time" in str(d).lower()), None)

    lat_name = find_dimension(dataset, ("lat", "latitude", "Latitude", "nlat"))
    lon_name = find_dimension(dataset, ("lon", "longitude", "Longitude", "nlon"))

    # Curvilinear grids (for example NSIDC sea ice) carry 2-D latitude and
    # longitude variables over y/x dimensions instead of lat/lon dimensions.
    lat2d_name = None
    lon2d_name = None
    if lat_name is None:
        for candidate in ("lat", "latitude", "Latitude"):
            if candidate in dataset.variables and dataset[candidate].ndim == 2:
                lat2d_name = candidate
                break
    if lon_name is None:
        for candidate in ("lon", "longitude", "Longitude"):
            if candidate in dataset.variables and dataset[candidate].ndim == 2:
                lon2d_name = candidate
                break

    notes = list(extra_notes or [])
    rows = []

    if mode == "grid":
        if lat_name is None or lon_name is None:
            fail("grid mode needs 1-D latitude and longitude dimension coordinates")
        latitude = dataset[lat_name].values
        longitude = dataset[lon_name].values
        if field.dims != (time_name, lat_name, lon_name) and field.dims != (lat_name, lon_name):
            fail(f"unexpected variable dimensions {field.dims} for grid mode")
        last = field.isel({time_name: -1}) if time_name in field.dims else field
        values = numpy.asarray(last.values, dtype=float)
        lat = numpy.asarray(latitude, dtype=float)
        lon = numpy.asarray(longitude, dtype=float)

        def edges(v):
            mid = 0.5 * (v[:-1] + v[1:])
            return numpy.concatenate(([v[0] - (mid[0] - v[0])], mid,
                                      [v[-1] + (v[-1] - mid[-1])]))

        lat_edges = edges(lat)
        lon_edges = edges(lon)
        for i in range(values.shape[0]):
            for j in range(values.shape[1]):
                value = float(values[i, j])
                if not numpy.isfinite(value):
                    continue
                rows.append((lat_edges[i], lat_edges[i + 1], lon_edges[j], lon_edges[j + 1], value))
        header = "lat_min,lat_max,lon_min,lon_max,concentration"
    else:
        if mode == "extent":
            if lat_name is not None and lon_name is not None:
                areas = cell_areas(dataset[lat_name].values, dataset[lon_name].values,
                                   lat_name, lon_name)
                area_dims = [lat_name, lon_name]
            elif lat2d_name is not None and lon2d_name is not None:
                areas = cell_areas_2d(dataset[lat2d_name].values, dataset[lon2d_name].values,
                                      dataset[lat2d_name].dims)
                area_dims = list(dataset[lat2d_name].dims)
                notes.append("curvilinear 2-D latitude/longitude grid")
            else:
                fail("extent mode needs latitude and longitude coordinates")
            # Sum cell areas where the concentration reaches the 15% threshold.
            series = (field >= threshold) * areas
            series = series.sum(dim=area_dims, skipna=True)
            units = units or "million km2"
            scale_output = 1.0 / 1e6
            notes.append(f"extent uses the {threshold:.2f} concentration threshold")
        else:
            if lat_name is not None and lon_name is not None and lat_name in field.dims and lon_name in field.dims:
                weights = numpy.cos(numpy.deg2rad(dataset[lat_name].values))
                weight_array = xarray.DataArray(weights, dims=(lat_name,))
                series = field.weighted(weight_array).mean(dim=[lat_name, lon_name], skipna=True)
                notes.append("cosine-latitude area weighting applied")
            elif lat2d_name is not None:
                weights = numpy.cos(numpy.deg2rad(dataset[lat2d_name].values))
                weight_array = xarray.DataArray(weights, dims=dataset[lat2d_name].dims)
                series = field.weighted(weight_array).mean(dim=list(weight_array.dims), skipna=True)
                notes.append("cosine-latitude weighting on a curvilinear grid")
            else:
                other = [d for d in field.dims if d != time_name]
                series = field.mean(dim=other, skipna=True) if other else field
            scale_output = 1.0

        if time_name is not None and time_name in series.dims:
            times = dataset[time_name].values
            try:
                import pandas as pd
                index = pd.DatetimeIndex(pd.to_datetime(times))
                if index.has_duplicates:
                    series = series.resample({time_name: "1MS"}).mean()
                    notes.append("duplicate timestamps averaged to monthly means")
                    times = series[time_name].values
            except Exception:
                pass
            years = decimal_years(times)
        else:
            years = numpy.arange(len(series))
            notes.append("no time coordinate found; index written instead of a year")

        values = numpy.asarray(series.values, dtype=float)
        count = min(len(years), len(values))
        for i in range(count):
            value = float(values[i]) * scale_output
            if not numpy.isfinite(value):
                continue
            rows.append((float(years[i]), value))
        header = "time,value"

    with open(out_csv, "w", encoding="utf-8") as handle:
        handle.write(header + "\n")
        for row in rows:
            handle.write(",".join(f"{value:.8g}" for value in row) + "\n")

    source_units = ""
    try:
        source_units = str(dataset[name].attrs.get("units", ""))
    except Exception:
        pass

    manifest = {
        "wrapper_version": WRAPPER_VERSION,
        "files": list(files),
        "variable": name,
        "source_units": source_units,
        "units": units or source_units,
        "mode": mode,
        "rows": len(rows),
        "output": out_csv,
        "notes": notes,
    }
    with open(out_csv + ".manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    print(json.dumps(manifest, indent=2))
    return manifest


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def command_search(args) -> int:
    granules = perform_search(args)
    summaries = [granule_summary(granule) for granule in granules]
    print(json.dumps({"count": len(summaries), "granules": summaries}, indent=2))
    return 0


def command_fetch(args) -> int:
    if args.dry_run:
        plan = {
            "action": "fetch",
            "short_name": args.short_name,
            "version": args.version,
            "start": args.start,
            "end": args.end,
            "bbox": args.bbox,
            "max_granules": args.max,
            "download_dir": args.out,
            "output_csv": args.out_csv,
            "variable": args.variable,
            "mode": args.mode,
        }
        print(json.dumps(plan, indent=2))
        return 0

    earthaccess = lazy_import("earthaccess", "Install earthaccess first.")
    login()
    granules = perform_search(args)
    print(f"matched {len(granules)} granule(s)", file=sys.stderr)
    if not granules:
        fail("no granules matched the search")
    os.makedirs(args.out, exist_ok=True)
    try:
        files = earthaccess.download(granules, args.out)
    except Exception as error:
        fail(f"download failed: {error}")
    if isinstance(files, str):
        files = [files]
    files = [path for path in files if path]
    if not files:
        fail("download returned no files")
    print(f"downloaded {len(files)} file(s)", file=sys.stderr)
    convert_files(files, args.out_csv, args.variable, args.mode, args.units,
                  args.threshold, extra_notes=["downloaded by earthaccess"])
    return 0


def command_convert(args) -> int:
    files = list(args.files)
    if not files:
        fail("provide one or more NetCDF files to convert")
    convert_files(files, args.out, args.variable, args.mode, args.units, args.threshold)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Search, download and convert NASA Earthdata granules to CSV."
    )
    parser.add_argument("--version-only", action="store_true", help="print the wrapper version")
    sub = parser.add_subparsers(dest="command")

    def add_search_arguments(target):
        target.add_argument("--short-name", required=True, help="CMR collection short name")
        target.add_argument("--version", default="", help="CMR collection version")
        target.add_argument("--concept-id", default="", help="CMR collection concept id")
        target.add_argument("--start", default="", help="start date YYYY-MM-DD")
        target.add_argument("--end", default="", help="end date YYYY-MM-DD")
        target.add_argument("--bbox", default="", help="west,south,east,north")
        target.add_argument("--max", type=int, default=0, help="limit the number of granules")

    search = sub.add_parser("search", help="search CMR and summarise granules")
    add_search_arguments(search)
    search.set_defaults(function=command_search)

    fetch = sub.add_parser("fetch", help="download granules and convert them to CSV")
    add_search_arguments(fetch)
    fetch.add_argument("--out", default="data/nasa", help="download directory")
    fetch.add_argument("--out-csv", required=True, help="CSV file to write")
    fetch.add_argument("--variable", default="", help="comma-separated variable candidates")
    fetch.add_argument("--mode", default="timeseries", choices=["timeseries", "extent", "grid"])
    fetch.add_argument("--units", default="", help="override the output units label")
    fetch.add_argument("--threshold", type=float, default=0.15)
    fetch.add_argument("--dry-run", action="store_true", help="print the plan without networking")
    fetch.set_defaults(function=command_fetch)

    convert = sub.add_parser("convert", help="convert local NetCDF files to CSV")
    convert.add_argument("files", nargs="+")
    convert.add_argument("--out", required=True, help="CSV file to write")
    convert.add_argument("--variable", default="", help="comma-separated variable candidates")
    convert.add_argument("--mode", default="timeseries", choices=["timeseries", "extent", "grid"])
    convert.add_argument("--units", default="", help="override the output units label")
    convert.add_argument("--threshold", type=float, default=0.15)
    convert.set_defaults(function=command_convert)

    return parser


def main(argv) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.version_only:
        print(WRAPPER_VERSION)
        return 0
    if not hasattr(args, "function"):
        parser.print_help()
        return 2
    return args.function(args)


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except KeyboardInterrupt:
        sys.exit(130)
