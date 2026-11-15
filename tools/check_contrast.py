#!/usr/bin/env python3
"""WCAG contrast audit for the workbench palette.

Reads the :root token block from src/ui_page.cpp, then reports the contrast of
every text token against the darkest and mid-dark surfaces the workbench uses.
Text tokens must reach 7:1 (WCAG AAA for body text); graphic tokens (series
lines, accents used as fills) must reach 3:1. Exits non-zero on a regression.
"""
import re
import sys

SOURCE = "src/ui_page.cpp"
SURFACES = {"--bg": "#0b0806", "--panel": "#1b130d", "--panel2": "#241a12",
            "sidebar": "#1a120c", "lab": "#0d0906", "card": "#120d08"}
TEXT_TOKENS = ["--text", "--muted"]
GRAPHIC_TOKENS = ["--teal", "--sky", "--amber", "--rose", "--green", "--violet",
                  "--series-obs", "--series-trend", "--series-season",
                  "--series-rem", "--series-alt", "--series-violet"]
TEXT_MIN = 7.0
GRAPHIC_MIN = 3.0


def luminance(hex_color):
    h = hex_color.lstrip("#")
    channels = [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    linear = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4 for c in channels]
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]


def contrast(a, b):
    la, lb = sorted([luminance(a), luminance(b)], reverse=True)
    return (la + 0.05) / (lb + 0.05)


def main():
    text = open(SOURCE, encoding="utf-8").read()
    block = re.search(r":root\{(.*?)\}", text, re.S)
    if not block:
        print("no :root block found", file=sys.stderr)
        return 1
    tokens = dict(re.findall(r"(--[\w-]+)\s*:\s*(#[0-9a-fA-F]{6})", block.group(1)))
    missing = [t for t in TEXT_TOKENS + GRAPHIC_TOKENS if t not in tokens]
    if missing:
        print("missing tokens: " + ", ".join(missing), file=sys.stderr)
        return 1

    failures = []
    print("%-16s %-9s %-9s %6s  %s" % ("token", "value", "surface", "ratio", "verdict"))
    for token in TEXT_TOKENS + GRAPHIC_TOKENS:
        floor = TEXT_MIN if token in TEXT_TOKENS else GRAPHIC_MIN
        for name, surface in SURFACES.items():
            ratio = contrast(tokens[token], surface)
            ok = ratio >= floor
            if not ok:
                failures.append((token, name, ratio, floor))
            print("%-16s %-9s %-9s %6.2f  %s" % (token, tokens[token], name, ratio,
                                                 "pass" if ok else "FAIL (needs %.1f)" % floor))
    # Bezel text: dark ink on the light bezel.
    bezel = tokens.get("--bezel", "#c8b294")
    ink = "#160f0a"
    ratio = contrast(ink, bezel)
    print("%-16s %-9s %-9s %6.2f  %s" % ("bezel ink", ink, bezel, ratio,
                                         "pass" if ratio >= TEXT_MIN else "FAIL"))
    if ratio < TEXT_MIN:
        failures.append(("bezel ink", "bezel", ratio, TEXT_MIN))

    if failures:
        print("\n%d contrast failure(s)" % len(failures), file=sys.stderr)
        return 1
    print("\nall palette contrasts pass (text >= %.1f:1, graphics >= %.1f:1)"
          % (TEXT_MIN, GRAPHIC_MIN))
    return 0


if __name__ == "__main__":
    sys.exit(main())
