#!/usr/bin/env bash
# Installs a Linux-hosted mingw-w64 cross-compiler without root, into .tools/mingw.
#
# The Windows binary is built with:
#   .tools/mingw/usr/bin/x86_64-w64-mingw32-g++-posix
# so nothing has to be installed system-wide and the setup is repeatable.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PREFIX="$ROOT/.tools/mingw"
WORK="$ROOT/.tools/debs"
SUITE="${SUITE:-noble}"

if [ -x "$PREFIX/usr/bin/x86_64-w64-mingw32-g++-posix" ]; then
    echo "toolchain already present at $PREFIX"
    exit 0
fi

mkdir -p "$WORK" "$PREFIX"
cd "$WORK"

for component in main universe; do
    if [ ! -f "$component.gz" ]; then
        curl -sSL -o "$component.gz" \
            "http://archive.ubuntu.com/ubuntu/dists/$SUITE/$component/binary-amd64/Packages.gz"
    fi
    gunzip -kf "$component.gz"
done

want="binutils-mingw-w64-x86-64 gcc-mingw-w64-base gcc-mingw-w64-x86-64 \
gcc-mingw-w64-x86-64-posix g++-mingw-w64-x86-64 g++-mingw-w64-x86-64-posix \
mingw-w64-common mingw-w64-x86-64-dev"

awk -v want="$want" '
BEGIN { n = split(want, names, " "); for (i = 1; i <= n; i++) targets[names[i]] = 1 }
/^Package: / { pkg = $2; keep = (pkg in targets) }
keep && /^Filename: / { print $2 }
' main universe | sort -u > files.txt

while read -r path; do
    name="$(basename "$path")"
    [ -f "$name" ] || curl -sSL -O "http://archive.ubuntu.com/ubuntu/$path"
    dpkg -x "$name" "$PREFIX"
done < files.txt

echo "installed: $("$PREFIX/usr/bin/x86_64-w64-mingw32-g++-posix" --version | head -1)"
