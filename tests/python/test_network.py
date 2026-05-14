"""
Test haplotype network algorithms.
"""

import pytest
from haplokit.network import (
    compute_tcs_network,
    compute_msn_network,
    compute_mjn_network,
    find_network_backend
)
from haplokit.plot import plot_hap_network


def test_find_backend():
    """Test backend discovery."""
    backend = find_network_backend()
    # May be None if not built yet
    if backend:
        assert backend.exists()


def test_msn_network():
    """Test MSN algorithm."""
    sequences = [
        "ATCGATCG",
        "ATCGTTCG",
        "TTCGATCG",
        "ATCGATCC"
    ]
    samples = [
        ["S1", "S2"],
        ["S3"],
        ["S4"],
        ["S5"]
    ]

    try:
        network = compute_msn_network(sequences, samples, epsilon=0)

        # Check structure
        assert "nodes" in network
        assert "edges" in network
        assert len(network["nodes"]) == 4

        # MSN should have n-1 edges (minimum spanning tree)
        assert len(network["edges"]) >= 3

        # Check node structure
        for node in network["nodes"]:
            assert "id" in node
            assert "sequence" in node
            assert "samples" in node
            assert "is_median" in node

        # Check edge structure
        for edge in network["edges"]:
            assert "source" in edge
            assert "target" in edge
            assert "weight" in edge

    except RuntimeError as e:
        if "not found" in str(e):
            pytest.skip("C++ backend not built")
        raise


def test_tcs_network():
    """Test TCS algorithm."""
    sequences = [
        "ATCG",
        "ATCC",
        "TTCG"
    ]
    samples = [
        ["S1"],
        ["S2"],
        ["S3"]
    ]

    try:
        network = compute_tcs_network(sequences, samples)

        assert "nodes" in network
        assert "edges" in network
        assert len(network["nodes"]) >= 3

    except RuntimeError as e:
        if "not found" in str(e):
            pytest.skip("C++ backend not built")
        raise


def test_mjn_network():
    """Test MJN algorithm."""
    sequences = [
        "ATCGATCG",
        "ATCGTTCG",
        "TTCGTTCG"
    ]
    samples = [
        ["S1"],
        ["S2"],
        ["S3"]
    ]

    try:
        network = compute_mjn_network(sequences, samples, epsilon=0)

        assert "nodes" in network
        assert "edges" in network

        # MJN may add median vertices
        assert len(network["nodes"]) >= 3

        # Check for median vertices
        has_median = any(node["is_median"] for node in network["nodes"])
        # May or may not have medians depending on topology

    except RuntimeError as e:
        if "not found" in str(e):
            pytest.skip("C++ backend not built")
        raise


def test_plot_hap_network_uses_cpp_network_result(tmp_path, monkeypatch):
    calls = {}

    def fake_compute(sequences, samples):
        calls["sequences"] = sequences
        calls["samples"] = samples
        return {
            "nodes": [
                {"id": 0, "sequence": sequences[0], "samples": ["S1", "S2"], "is_median": False},
                {"id": 1, "sequence": sequences[1], "samples": ["S3"], "is_median": False},
                {"id": 2, "sequence": "MEDIAN", "samples": [], "is_median": True},
            ],
            "edges": [
                {"source": 0, "target": 2, "weight": 1},
                {"source": 2, "target": 1, "weight": 1},
            ],
        }

    monkeypatch.setattr("haplokit._network.compute_tcs_network", fake_compute)

    output = plot_hap_network(
        ["H001", "H002"],
        ["ATCG", "ATCC"],
        [2, 1],
        tmp_path / "network.png",
        pop_data={"H001": [("Pop1", 2)], "H002": [("Pop2", 1)]},
    )

    assert calls["sequences"] == ["ATCG", "ATCC"]
    assert calls["samples"] == [["H001", "H001"], ["H002"]]
    assert output.exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
