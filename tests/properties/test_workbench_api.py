"""Hypothesis property tests for the workbench HTTP API.

The server is started on demand (or reused if it already listens) and every
property drives the real endpoints: no mocks, no recorded fixtures.

    .venv-test/bin/pytest tests/properties -q
"""

import csv
import io
import json
import os
import socket
import subprocess
import time
import urllib.error
import urllib.request

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PORT = int(os.environ.get("CLIMATE_TEST_PORT", "3399"))
BASE = "http://127.0.0.1:%d" % PORT


def _listening():
    with socket.socket() as sock:
        sock.settimeout(0.3)
        return sock.connect_ex(("127.0.0.1", PORT)) == 0


@pytest.fixture(scope="session", autouse=True)
def server():
    process = None
    if not _listening():
        process = subprocess.Popen([os.path.join(ROOT, "climate"), "ui", "--port", str(PORT),
                                    "--host", "127.0.0.1"], cwd=ROOT,
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for _ in range(120):
            if _listening():
                break
            time.sleep(0.1)
        else:
            process.kill()
            pytest.fail("the workbench did not start on port %d" % PORT)
    yield
    if process is not None:
        process.terminate()
        process.wait(timeout=10)


def get(path):
    with urllib.request.urlopen(BASE + path, timeout=60) as response:
        return response.status, json.loads(response.read().decode())


def post(path, body):
    request = urllib.request.Request(BASE + path, data=json.dumps(body).encode(),
                                     headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            return response.status, json.loads(response.read().decode())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read().decode() or "{}")


@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(seed=st.integers(min_value=0, max_value=10_000), count=st.integers(min_value=1, max_value=6))
def test_batch_returns_one_result_per_formula(seed, count):
    status, catalog = get("/api/lab")
    assert status == 200
    formulas = [f["id"] for group in catalog["groups"] for f in group["formulas"]]
    chosen = [formulas[(seed * 37 + i * 101) % len(formulas)] for i in range(count)]
    status, body = post("/api/lab/run", {"calls": [{"id": name} for name in chosen]})
    assert status == 200
    assert body["requested"] == count
    assert len(body["results"]) == count
    for result in body["results"]:
        assert isinstance(result["ok"], bool)
        assert ("result" in result) or ("error" in result)


@settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(slope=st.floats(min_value=-0.4, max_value=0.4, allow_nan=False),
       noise=st.floats(min_value=0.0, max_value=2.0, allow_nan=False),
       n=st.integers(min_value=24, max_value=90))
def test_analysis_never_fails_with_a_server_error(slope, noise, n):
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(["time", "value"])
    for i in range(n):
        writer.writerow(["%.4f" % (2000 + i / 12.0), 14.0 + slope * i + noise * ((i * 7) % 5)])
    status, body = post("/api/analyze", {"csv": buffer.getvalue(), "period": 12, "alpha": 0.05})
    assert status in (200, 400), body
    if status == 200:
        assert body["samples"] >= 3
        assert "trend" in body and "ensemble" in body["trend"]


@settings(max_examples=10, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(lat=st.floats(min_value=-90, max_value=90, allow_nan=False))
def test_power_pull_validates_latitude_without_crashing(lat):
    status, body = post("/api/live/pull", {"source": "power", "lat": lat, "lon": 10,
                                           "start": "2010-01-01", "end": "2010-12-31",
                                           "parameter": "T2M"})
    assert status == 200
    assert isinstance(body.get("ok"), bool)
    if not body["ok"]:
        assert body["error"]


def test_live_status_lists_public_endpoints_only():
    status, body = get("/api/live/status")
    assert status == 200
    assert body["sources"]
    for source in body["sources"]:
        assert source["url"].startswith("https://")
        assert "password" not in json.dumps(source).lower()
