"""
Python prototype for network algorithms (fallback when C++ backend unavailable).
"""

import numpy as np
from typing import List, Dict, Any, Tuple
from collections import defaultdict
import heapq


def hamming_distance(seq1: str, seq2: str) -> int:
    """Compute Hamming distance between two sequences."""
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def compute_msn_python(
    sequences: List[str],
    samples: List[List[str]],
    epsilon: int = 0
) -> Dict[str, Any]:
    """
    Python implementation of MSN (Minimum Spanning Network).
    Uses Kruskal's algorithm with union-find.
    """
    n = len(sequences)

    # Union-Find
    parent = list(range(n))
    rank = [0] * n

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        root_x, root_y = find(x), find(y)
        if root_x == root_y:
            return False
        if rank[root_x] < rank[root_y]:
            parent[root_x] = root_y
        elif rank[root_x] > rank[root_y]:
            parent[root_y] = root_x
        else:
            parent[root_y] = root_x
            rank[root_x] += 1
        return True

    # Compute all pairwise distances
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            dist = hamming_distance(sequences[i], sequences[j])
            edges.append((dist, i, j))

    # Sort by distance
    edges.sort()

    # Kruskal's algorithm
    result_edges = []
    components = n
    max_threshold = float('inf')

    for dist, u, v in edges:
        if dist > max_threshold:
            break

        if union(u, v):
            result_edges.append((u, v, dist))
            components -= 1

            if components == 1 and max_threshold == float('inf'):
                max_threshold = dist + epsilon

    # Build output
    nodes = [
        {
            "id": i,
            "sequence": sequences[i],
            "samples": samples[i] if i < len(samples) else [],
            "is_median": False
        }
        for i in range(n)
    ]

    edges_out = [
        {"source": u, "target": v, "weight": w}
        for u, v, w in result_edges
    ]

    return {"nodes": nodes, "edges": edges_out}


def compute_tcs_python(
    sequences: List[str],
    samples: List[List[str]]
) -> Dict[str, Any]:
    """
    Python implementation of TCS (simplified).
    Connects haplotypes at increasing distance thresholds.
    """
    n = len(sequences)

    # Group pairs by distance
    pairs_by_dist = defaultdict(list)
    for i in range(n):
        for j in range(i + 1, n):
            dist = hamming_distance(sequences[i], sequences[j])
            pairs_by_dist[dist].append((i, j))

    # Track components
    component_id = list(range(n))
    edges = []

    # Process pairs in order of increasing distance
    for dist in sorted(pairs_by_dist.keys()):
        comp_a, comp_b = -1, -1

        for u, v in pairs_by_dist[dist]:
            comp_u, comp_v = component_id[u], component_id[v]

            if comp_u == comp_v:
                continue

            if comp_u > comp_v:
                comp_u, comp_v = comp_v, comp_u
                u, v = v, u

            if comp_a < 0:
                comp_a, comp_b = comp_u, comp_v

            if comp_u == comp_a and comp_v == comp_b:
                edges.append((u, v, dist))

        # Merge components
        if comp_a >= 0:
            for i in range(n):
                if component_id[i] == comp_b:
                    component_id[i] = comp_a
                elif component_id[i] > comp_b:
                    component_id[i] -= 1

    # Build output
    nodes = [
        {
            "id": i,
            "sequence": sequences[i],
            "samples": samples[i] if i < len(samples) else [],
            "is_median": False
        }
        for i in range(n)
    ]

    edges_out = [
        {"source": u, "target": v, "weight": w}
        for u, v, w in edges
    ]

    return {"nodes": nodes, "edges": edges_out}


def compute_mjn_python(
    sequences: List[str],
    samples: List[List[str]],
    epsilon: int = 0
) -> Dict[str, Any]:
    """
    Python implementation of MJN (simplified).
    Starts with MSN, then adds median vertices.
    """
    # Start with MSN
    return compute_msn_python(sequences, samples, epsilon)
