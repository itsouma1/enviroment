"""BDD steps for tests/features/change_is_always_occurring.feature.

The scenarios describe the product premise in plain language and drive the same
HTTP endpoints the browser uses, so a change in behaviour shows up as a failing
sentence rather than a failing assertion deep in a helper.

    .venv-test/bin/pytest tests/bdd -q
"""

import io
import json
import os
import socket
import subprocess
import time
import urllib.error
import urllib.request

import pytest
from pytest_bdd import given, parsers, scenarios, then, when

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PORT = int(os.environ.get("CLIMATE_TEST_PORT", "3400"))
BASE = "http://127.0.0.1:%d" % PORT

scenarios("../features/change_is_always_occurring.feature")


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


def post(path, body):
    request = urllib.request.Request(BASE + path, data=json.dumps(body).encode(),
                                     headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=180) as response:
            return response.status, json.loads(response.read().decode())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read().decode() or "{}")


def get(path):
    with urllib.request.urlopen(BASE + path, timeout=180) as response:
        return response.status, json.loads(response.read().decode())


@pytest.fixture
def context():
    return {}


def _series(slope):
    buffer = io.StringIO()
    buffer.write("time,value\n")
    for i in range(180):
        buffer.write("%.4f,%.6f\n" % (2000 + i / 12.0, 14.0 + slope * i + 0.3 * ((i * 3) % 4)))
    return buffer.getvalue()


@given("a monthly series warming at 0.03 per month")
def warming_series(context):
    context["csv"] = _series(0.03)


@given("a monthly series with no trend")
def flat_series(context):
    context["csv"] = _series(0.0)


@when("the workbench analyses it")
def analyse(context):
    status, body = post("/api/analyze", {"csv": context["csv"], "period": 12, "alpha": 0.05})
    assert status == 200, body
    context["analysis"] = body


@then("a positive trend slope is reported")
def positive_slope(context):
    slope = context["analysis"]["trend"]["slope"]
    assert slope > 0.0, slope


@then("the trend is declared significant")
def significant(context):
    assert context["analysis"]["trend"]["significant"] is True


@then("the trend estimator ensemble reports a primary method")
def ensemble(context):
    ensemble = context["analysis"]["trend"]["ensemble"]
    assert ensemble["primary"]["name"]
    assert ensemble["primary"]["slope"] is not None


@then("the trend is not declared significant")
def not_significant(context):
    assert context["analysis"]["trend"]["significant"] is False


@when("the spatial grid is requested")
def spatial(context):
    status, body = get("/api/grid-trend")
    assert status == 200
    context["grid"] = body


@then("every cell has a trend slope")
def cells(context):
    grid = context["grid"]
    assert len(grid["slopes"]) == grid["rows"] * grid["columns"]


@then("at least one region is flagged as a hotspot")
def hotspots(context):
    assert context["grid"]["hotspot_count"] >= 1


@given("the formula lab catalogue")
def catalogue(context):
    status, body = get("/api/lab")
    assert status == 200
    context["lab"] = body


@when("a batch of three formulas is run")
def batch(context):
    formulas = [f["id"] for group in context["lab"]["groups"] for f in group["formulas"]]
    chosen = [name for name in ("physics::evi", "trend::sens_slope", "spatial::getis_ord_gi")
              if name in formulas] or formulas[:3]
    status, body = post("/api/lab/run", {"calls": [{"id": name} for name in chosen]})
    assert status == 200
    context["batch"] = body


@then("three results are returned")
def three_results(context):
    assert context["batch"]["requested"] == len(context["batch"]["results"]) == 3


@then("each result either has a value or a documented error")
def result_shape(context):
    for result in context["batch"]["results"]:
        assert result.get("ok") in (True, False)
        assert ("result" in result) or ("error" in result)


@when("the GISTEMP series is pulled")
def pull_gistemp(context):
    status, body = post("/api/live/pull", {"source": "gistemp"})
    assert status == 200
    context["live"] = body


@then("more than a thousand samples are returned")
def many_samples(context):
    live = context["live"]
    if not live.get("ok"):
        pytest.skip("network unavailable: %s" % live.get("error"))
    assert live["series"]["samples"] > 1000


@then("the provenance records the NASA endpoint")
def provenance(context):
    live = context["live"]
    if not live.get("ok"):
        pytest.skip("network unavailable")
    assert "giss.nasa.gov" in live["series"]["url"]
    assert live["series"]["retrieved_at"]
