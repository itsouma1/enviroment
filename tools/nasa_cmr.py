#!/usr/bin/env python3
"""Query the NASA CMR API for granule metadata without downloading data.

This uses only the Python standard library. It needs no Earthdata Login account:
CMR search results are public metadata. Data access still needs credentials.
"""

import argparse
import json
import urllib.error
import urllib.parse
import urllib.request

CMR_GRANULES = "https://cmr.earthdata.nasa.gov/search/granules.json"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--short-name", required=True)
    parser.add_argument("--start", default="")
    parser.add_argument("--end", default="")
    parser.add_argument("--max", type=int, default=5)
    parser.add_argument("--version", default="")
    args = parser.parse_args()

    params = {
        "short_name": args.short_name,
        "page_size": max(1, min(args.max, 50)),
        "sort_key[]": "-start_date",
    }
    if args.start and args.end:
        params["temporal"] = f"{args.start},{args.end}"
    elif args.start:
        params["temporal"] = f"{args.start},"
    if args.version:
        params["version"] = args.version

    url = CMR_GRANULES + "?" + urllib.parse.urlencode(params, doseq=True)
    request = urllib.request.Request(
        url,
        headers={"Accept": "application/json", "User-Agent": "climate-workbench/1.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=25) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except (urllib.error.URLError, TimeoutError, ValueError) as exc:
        print(json.dumps({"ok": False, "url": url, "error": str(exc)}))
        return 2

    feed = payload.get("feed", {})
    entries = feed.get("entry", [])
    if isinstance(entries, dict):
        entries = [entries]
    items = []
    for entry in entries[: args.max]:
        links = []
        for link in entry.get("links", []):
            href = link.get("href", "")
            if href and ("opendap" in href or href.endswith(".nc") or "data" in href):
                links.append(href)
        items.append(
            {
                "title": entry.get("title", ""),
                "id": entry.get("id", ""),
                "time_start": entry.get("time_start", ""),
                "time_end": entry.get("time_end", ""),
                "links": links[:3],
            }
        )

    print(
        json.dumps(
            {
                "ok": True,
                "url": url,
                "hits": feed.get("hits", len(items)),
                "items": items,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
