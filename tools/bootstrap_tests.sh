#!/usr/bin/env bash
# Install everything the extended test suites need, reproducibly.
#
#   tools/bootstrap_tests.sh
#
# Creates .venv-test (cmake, ctest, gcovr, pytest, pytest-bdd, hypothesis,
# earthaccess, xarray, netCDF4) and vendors GoogleTest and RapidCheck under
# .tools/deps. Nothing here is required to build or run the workbench itself.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

python="${PYTHON:-python3}"

if [ ! -d .venv-test ]; then
  "$python" -m venv .venv-test
fi
.venv-test/bin/pip install --quiet --upgrade pip
.venv-test/bin/pip install --quiet cmake gcovr hypothesis pytest pytest-bdd \
  earthaccess xarray netCDF4 pandas requests

mkdir -p .tools/deps
if [ ! -d .tools/deps/googletest ]; then
  git clone --depth 1 https://github.com/google/googletest.git .tools/deps/googletest
fi
if [ ! -d .tools/deps/rapidcheck ]; then
  git clone --depth 1 https://github.com/emil-e/rapidcheck.git .tools/deps/rapidcheck
fi

echo "cmake:      $(.venv-test/bin/cmake --version | head -1)"
echo "gtest:      $(ls .tools/deps/googletest | head -1)"
echo "rapidcheck: $(ls .tools/deps/rapidcheck | head -1)"
echo "pytest:     $(.venv-test/bin/pytest --version)"
echo
echo "Next: make tests-cmake && make coverage && make property && make bdd"
