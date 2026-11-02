"""Unit tests for src/nvidia_client.py - all HTTP is mocked, nothing goes out.

The fake session below mimics the behaviours we actually observed from the
hosted NIMs on 2026-11-02: 200 with `Nvcf-Status: fulfilled`, 422 with a
one-line `{"error": ...}` validation message, 429/5xx, and the 202 +
`nvcf-reqid` polling pattern.
"""
from __future__ import annotations

import json

import pytest
import requests

import nvidia_client
from nvidia_client import (
    NvidiaApiError,
    NvidiaAuthError,
    NvidiaClient,
    parse_mpnn_mfasta,
)


class FakeResponse:
    def __init__(self, status_code=200, payload=None, headers=None, text=None):
        self.status_code = status_code
        self._payload = payload
        self.headers = headers or {}
        self._text = text

    @property
    def text(self):
        if self._text is not None:
            return self._text
        return json.dumps(self._payload)

    def json(self):
        if self._payload is None:
            raise ValueError("no json")
        return self._payload


class FakeSession:
    """Returns queued responses and records the requests it was given."""

    def __init__(self, responses):
        self.responses = list(responses)
        self.posts = []
        self.gets = []

    def post(self, url, headers=None, json=None, timeout=None):
        self.posts.append({"url": url, "headers": headers, "json": json,
                           "timeout": timeout})
        return self._next()

    def get(self, url, headers=None, timeout=None):
        self.gets.append({"url": url, "headers": headers, "timeout": timeout})
        return self._next()

    def _next(self):
        item = self.responses.pop(0)
        if isinstance(item, Exception):
            raise item
        return item


@pytest.fixture(autouse=True)
def _api_key(monkeypatch):
    monkeypatch.setenv("NVIDIA_API_KEY", "test-key-not-a-real-key")


def make_client(responses, **kwargs):
    session = FakeSession(responses)
    client = NvidiaClient(
        session=session, sleep=lambda _s: None, max_retries=kwargs.pop("max_retries", 2),
        **kwargs,
    )
    return client, session


# --- authentication --------------------------------------------------------

def test_missing_key_raises_auth_error(monkeypatch):
    monkeypatch.delenv("NVIDIA_API_KEY", raising=False)
    with pytest.raises(NvidiaAuthError) as excinfo:
        NvidiaClient()
    assert "NVIDIA_API_KEY" in str(excinfo.value)


def test_blank_key_raises_auth_error(monkeypatch):
    monkeypatch.setenv("NVIDIA_API_KEY", "   ")
    with pytest.raises(NvidiaAuthError):
        NvidiaClient()


def test_key_is_sent_as_bearer_and_not_in_the_body():
    client, session = make_client([FakeResponse(payload={"ok": True})])
    client.post("/v1/thing", {"a": 1})
    headers = session.posts[0]["headers"]
    assert headers["Authorization"] == "Bearer test-key-not-a-real-key"
    assert "test-key" not in json.dumps(session.posts[0]["json"])


def test_401_is_an_auth_error_without_echoing_the_key():
    client, _ = make_client(
        [FakeResponse(401, {"error": "invalid token"})], max_retries=0
    )
    with pytest.raises(NvidiaAuthError) as excinfo:
        client.post("/v1/thing", {})
    message = str(excinfo.value)
    assert "401" in message and "invalid token" in message
    assert "test-key-not-a-real-key" not in message


# --- retries and polling ---------------------------------------------------

def test_429_is_retried_then_succeeds():
    client, session = make_client(
        [
            FakeResponse(429, {"error": "rate limit"}, {"Retry-After": "0"}),
            FakeResponse(payload={"output_pdb": "ATOM"}),
        ]
    )
    assert client.post("/v1/thing", {}) == {"output_pdb": "ATOM"}
    assert len(session.posts) == 2


def test_5xx_is_retried_until_retries_run_out():
    client, session = make_client(
        [FakeResponse(503, {"error": "upstream"})] * 3, max_retries=2
    )
    with pytest.raises(NvidiaApiError) as excinfo:
        client.post("/v1/thing", {})
    assert excinfo.value.status_code == 503
    assert len(session.posts) == 3


def test_422_is_not_retried_and_surfaces_the_validation_message():
    client, session = make_client(
        [FakeResponse(422, {"error": "Validation Errors: contigs: Field required"})]
    )
    with pytest.raises(NvidiaApiError) as excinfo:
        client.post("/v1/thing", {})
    assert "contigs: Field required" in str(excinfo.value)
    assert len(session.posts) == 1


def test_connection_errors_are_retried():
    client, session = make_client(
        [requests.ConnectionError("boom"), FakeResponse(payload={"ok": 1})]
    )
    assert client.post("/v1/thing", {}) == {"ok": 1}
    assert len(session.posts) == 2


def test_202_is_polled_until_fulfilled():
    client, session = make_client(
        [
            FakeResponse(202, None, {"nvcf-reqid": "abc-123"}, text=""),
            FakeResponse(202, None, {}, text=""),
            FakeResponse(payload={"output_pdb": "ATOM 1"}),
        ]
    )
    assert client.post("/v1/thing", {}) == {"output_pdb": "ATOM 1"}
    assert session.gets[0]["url"].endswith("/v2/nvcf/pexec/status/abc-123")
    assert len(session.gets) == 2


def test_202_without_a_request_id_is_an_error():
    client, _ = make_client([FakeResponse(202, None, {}, text="")])
    with pytest.raises(NvidiaApiError, match="no nvcf-reqid"):
        client.post("/v1/thing", {})


def test_polling_gives_up_after_the_poll_timeout():
    client, _ = make_client(
        [FakeResponse(202, None, {"nvcf-reqid": "abc"}, text="")]
        + [FakeResponse(202, None, {}, text="")] * 5,
        poll_timeout=0.0,
    )
    with pytest.raises(NvidiaApiError, match="still running"):
        client.post("/v1/thing", {})


def test_poll_error_status_is_raised():
    client, _ = make_client(
        [
            FakeResponse(202, None, {"nvcf-reqid": "abc"}, text=""),
            FakeResponse(500, {"error": "died"}),
        ]
    )
    with pytest.raises(NvidiaApiError, match="died"):
        client.post("/v1/thing", {})


def test_non_object_body_is_rejected():
    client, _ = make_client([FakeResponse(payload=["not", "an", "object"])])
    with pytest.raises(NvidiaApiError, match="expected a"):
        client.post("/v1/thing", {})


def test_backoff_prefers_retry_after_then_grows():
    client, _ = make_client([FakeResponse(payload={})])
    assert client._backoff(1, FakeResponse(429, {}, {"Retry-After": "7"})) == 7.0
    assert client._backoff(1, None) < client._backoff(4, None)


# --- RFdiffusion -----------------------------------------------------------

def test_rfdiffusion_request_body_and_result():
    client, session = make_client(
        [FakeResponse(payload={"output_pdb": "ATOM      1  N   GLY A   1",
                               "elapsed_ms": 2700})]
    )
    pdb = client.rfdiffusion_generate(
        input_pdb="ATOM  target",
        contigs="A1-180/0 C1-10/0 70-70",
        hotspot_res=["C4", "C5", "C6"],
        diffusion_steps=15,
        random_seed=1,
    )
    assert pdb.startswith("ATOM")
    body = session.posts[0]["json"]
    assert body == {
        "input_pdb": "ATOM  target",
        "contigs": "A1-180/0 C1-10/0 70-70",
        "diffusion_steps": 15,
        "hotspot_res": ["C4", "C5", "C6"],
        "random_seed": 1,
    }
    assert session.posts[0]["url"].endswith(nvidia_client.RFDIFFUSION_PATH)


def test_rfdiffusion_omits_optional_fields_when_unset():
    client, session = make_client([FakeResponse(payload={"output_pdb": "ATOM x"})])
    client.rfdiffusion_generate("pdb", "70-70")
    body = session.posts[0]["json"]
    assert "hotspot_res" not in body and "random_seed" not in body


def test_rfdiffusion_rejects_a_response_without_coordinates():
    client, _ = make_client([FakeResponse(payload={"output_pdb": ""})])
    with pytest.raises(NvidiaApiError, match="output_pdb"):
        client.rfdiffusion_generate("pdb", "70-70")


# --- ProteinMPNN -----------------------------------------------------------

MFASTA = (
    ">input, score=2.69, fixed_chains=['B'], designed_chains=['A']\n"
    "GGGGG\n"
    ">T=0.15, sample=1, score=1.0197, seq_recovery=0.1143\n"
    "STSSS\n"
    ">T=0.15, sample=2, score=1.0179, seq_recovery=0.1143\n"
    "SSSAS\n"
)


def test_proteinmpnn_request_body_and_parsing():
    client, session = make_client(
        [FakeResponse(payload={"mfasta": MFASTA, "scores": [1.02, 1.01],
                               "probs": [[[0.0]]]})]
    )
    result = client.proteinmpnn_design(
        input_pdb="ATOM backbone",
        input_pdb_chains=["A"],
        num_seq_per_target=2,
        sampling_temp=[0.15],
        random_seed=7,
    )
    body = session.posts[0]["json"]
    assert body["input_pdb_chains"] == ["A"]
    assert body["num_seq_per_target"] == 2
    assert body["sampling_temp"] == [0.15]
    assert body["use_soluble_model"] is True
    assert body["ca_only"] is False
    assert [s["sequence"] for s in result["sequences"]] == ["STSSS", "SSSAS"]
    assert result["sequences"][0]["score"] == pytest.approx(1.0197)
    assert result["sequences"][1]["sample"] == 2


def test_parse_mpnn_mfasta_drops_the_echoed_input_record():
    records = parse_mpnn_mfasta(MFASTA)
    assert len(records) == 2
    assert all("GGGGG" not in r["sequence"] for r in records)


def test_parse_mpnn_mfasta_tolerates_missing_score_fields():
    records = parse_mpnn_mfasta(">T=0.1, sample=x\nACDE\n")
    assert records[0]["sequence"] == "ACDE"
    assert records[0]["score"] is None
    assert records[0]["sample"] is None


def test_proteinmpnn_rejects_a_response_without_mfasta():
    client, _ = make_client([FakeResponse(payload={"scores": []})])
    with pytest.raises(NvidiaApiError, match="mfasta"):
        client.proteinmpnn_design("ATOM x")


# --- Boltz-2 ---------------------------------------------------------------

BOLTZ_OK = {
    "structures": [{"structure": "data_x", "format": "mmcif", "name": "n",
                    "source": "s"}],
    "iptm_scores": [0.79],
    "pair_chains_iptm_scores": [{"0": {"0": 0.96, "1": 0.1}, "1": {"0": 0.1, "1": 0.9}}],
}


def test_boltz2_request_body():
    client, session = make_client([FakeResponse(payload=BOLTZ_OK)])
    client.boltz2_predict(
        [
            {"id": "A", "molecule_type": "protein", "sequence": "acdef"},
            {"id": "D", "molecule_type": "protein", "sequence": "GHIKL"},
        ],
        recycling_steps=1,
        sampling_steps=30,
        without_potentials=True,
    )
    body = session.posts[0]["json"]
    assert body["polymers"][0] == {"molecule_type": "protein", "sequence": "ACDEF",
                                   "id": "A"}
    assert body["recycling_steps"] == 1
    assert body["sampling_steps"] == 30
    assert body["without_potentials"] is True
    assert session.posts[0]["url"].endswith(nvidia_client.BOLTZ2_PATH)


def test_boltz2_rejects_empty_and_oversized_input():
    client, _ = make_client([FakeResponse(payload=BOLTZ_OK)] * 3)
    with pytest.raises(ValueError, match="at least one polymer"):
        client.boltz2_predict([])
    with pytest.raises(ValueError, match="at most 12 polymers"):
        client.boltz2_predict(
            [{"molecule_type": "protein", "sequence": "A"}] * 13
        )
    with pytest.raises(ValueError, match="4096"):
        client.boltz2_predict(
            [{"id": "A", "molecule_type": "protein", "sequence": "A" * 4097}]
        )
    with pytest.raises(ValueError, match="non-empty"):
        client.boltz2_predict([{"id": "A", "molecule_type": "protein",
                                "sequence": "  "}])


def test_boltz2_rejects_a_response_without_structures():
    client, _ = make_client([FakeResponse(payload={"iptm_scores": [0.5]})])
    with pytest.raises(NvidiaApiError, match="no structures"):
        client.boltz2_predict([{"id": "A", "molecule_type": "protein",
                                "sequence": "ACDEF"}])


def test_module_level_helpers_use_one_client(monkeypatch):
    seen = {}

    class Stub:
        def __init__(self, **kwargs):
            seen["kwargs"] = kwargs

        def rfdiffusion_generate(self, input_pdb, contigs, **kw):
            seen["call"] = ("rfdiffusion", contigs, kw)
            return "ATOM"

    monkeypatch.setattr(nvidia_client, "NvidiaClient", Stub)
    assert nvidia_client.rfdiffusion_generate(
        "pdb", "70-70", timeout=5.0, diffusion_steps=20
    ) == "ATOM"
    assert seen["kwargs"] == {"timeout": 5.0}
    assert seen["call"][2] == {"diffusion_steps": 20}
