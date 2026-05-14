"""
Haplotype network algorithms interface.

Calls C++ backend for TCS, MSN, and MJN network computation.
"""

import json
import subprocess
from pathlib import Path
from typing import List, Dict, Any, Optional
import os


def find_network_backend() -> Optional[Path]:
    """Find haplokit_network_backend executable."""
    # Check environment variable
    if "HAPLOKIT_NETWORK_BIN" in os.environ:
        backend = Path(os.environ["HAPLOKIT_NETWORK_BIN"])
        if backend.exists():
            return backend

    # Check packaged binary
    package_dir = Path(__file__).parent
    backend = package_dir / "_bin" / "haplokit_network_backend"
    if backend.exists():
        return backend

    # Check build directory
    repo_root = package_dir.parent
    for build_dir in ["build-wsl", "build", "."]:
        backend = repo_root / build_dir / "haplokit_network_backend"
        if backend.exists():
            return backend
        # Windows executable
        backend_exe = repo_root / build_dir / "haplokit_network_backend.exe"
        if backend_exe.exists():
            return backend_exe

    return None


def compute_network(
    sequences: List[str],
    samples: List[List[str]],
    algorithm: str = "tcs",
    epsilon: int = 0,
    use_python_fallback: bool = True
) -> Dict[str, Any]:
    """
    Compute haplotype network using C++ backend or Python fallback.

    Args:
        sequences: List of haplotype sequences
        samples: List of sample lists for each haplotype
        algorithm: Network algorithm ("tcs", "msn", or "mjn")
        epsilon: Epsilon parameter for MSN/MJN relaxation
        use_python_fallback: Use Python implementation if C++ backend unavailable

    Returns:
        Dictionary with "nodes" and "edges" keys
    """
    backend = find_network_backend()

    # Try C++ backend first
    if backend:
        # Prepare input JSON
        input_data = {
            "sequences": sequences,
            "samples": samples,
            "epsilon": epsilon
        }
        input_json = json.dumps(input_data)

        # Call C++ backend
        try:
            # Check if backend is Linux ELF (needs WSL on Windows)
            import platform
            if platform.system() == "Windows":
                # Use WSL to run Linux binary
                result = subprocess.run(
                    ["wsl", str(backend).replace("\\", "/").replace("F:", "/mnt/f"), algorithm],
                    input=input_json,
                    capture_output=True,
                    text=True,
                    check=True
                )
            else:
                result = subprocess.run(
                    [str(backend), algorithm],
                    input=input_json,
                    capture_output=True,
                    text=True,
                    check=True
                )
            return json.loads(result.stdout)

        except (subprocess.CalledProcessError, json.JSONDecodeError, FileNotFoundError) as e:
            if not use_python_fallback:
                raise RuntimeError(f"Network computation failed: {e}") from e
            # Fall through to Python implementation

    # Python fallback
    if use_python_fallback:
        import sys
        from pathlib import Path
        archive_path = Path(__file__).parent.parent / "archive" / "python_reference_implementation"
        sys.path.insert(0, str(archive_path))
        from network_python import (
            compute_msn_python,
            compute_tcs_python,
            compute_mjn_python
        )

        if algorithm == "msn":
            return compute_msn_python(sequences, samples, epsilon)
        elif algorithm == "tcs":
            return compute_tcs_python(sequences, samples)
        elif algorithm == "mjn":
            return compute_mjn_python(sequences, samples, epsilon)
        else:
            raise ValueError(f"Unknown algorithm: {algorithm}")

    raise RuntimeError(
        "haplokit_network_backend not found and Python fallback disabled. "
        "Set HAPLOKIT_NETWORK_BIN or rebuild with CMake."
    )


def compute_tcs_network(
    sequences: List[str],
    samples: List[List[str]]
) -> Dict[str, Any]:
    """
    Compute TCS (Templeton-Crandall-Sing) statistical parsimony network.

    Args:
        sequences: List of haplotype sequences
        samples: List of sample lists for each haplotype

    Returns:
        Network with nodes and edges
    """
    return compute_network(sequences, samples, algorithm="tcs")


def compute_msn_network(
    sequences: List[str],
    samples: List[List[str]],
    epsilon: int = 0
) -> Dict[str, Any]:
    """
    Compute MSN (Minimum Spanning Network).

    Args:
        sequences: List of haplotype sequences
        samples: List of sample lists for each haplotype
        epsilon: Relaxation parameter (default 0 = strict MST)

    Returns:
        Network with nodes and edges
    """
    return compute_network(sequences, samples, algorithm="msn", epsilon=epsilon)


def compute_mjn_network(
    sequences: List[str],
    samples: List[List[str]],
    epsilon: int = 0
) -> Dict[str, Any]:
    """
    Compute MJN (Median-Joining Network).

    Iteratively adds median vertices to reduce network length.
    May create reticulation (loops) representing recombination or parallel mutation.

    Args:
        sequences: List of haplotype sequences
        samples: List of sample lists for each haplotype
        epsilon: Cost threshold for adding median vertices

    Returns:
        Network with nodes and edges (may include median vertices)
    """
    return compute_network(sequences, samples, algorithm="mjn", epsilon=epsilon)
