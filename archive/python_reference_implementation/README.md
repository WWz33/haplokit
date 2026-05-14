# Python Reference Implementation (Archive)

This directory contains archived **pure Python** haplotype network code and validation scripts.

## Status: ARCHIVED

⚠️ These files are **not production path**.

Production architecture now is:
- **C++ backend** computes haplotype networks
- **Python frontend** renders network figures

## Purpose

This archive is kept for:
- **Reference implementation** of MSN / TCS / MJN Python algorithms
- **C++ reproduction check** against earlier Python behavior
- **Visualization validation** with historical demo / e2e scripts
- **Educational reading** of earlier algorithm prototypes

## Files

- `network_python.py` - Pure Python implementations of `compute_msn_python()`, `compute_tcs_python()`, `compute_mjn_python()`
- `test_e2e_python_backend.py` - End-to-end test for Python algorithm backend + PopART-style plotting
- `demo_algorithms_popart.py` - Algorithm comparison demo for MSN / TCS / MJN with PopART-style figure output

## Notes

The Python algorithm scripts above are preserved because they were used as reference when reproducing behavior in C++ backend.

Package active path no longer uses these files directly.

## Production path

Current active path is:
- network compute: `haplokit/network.py`
- C++ backend: `src/cpp/network/`
- Python plotting frontend: `haplokit/_network.py`, `haplokit/plot.py`

## Usage (archive only)

```python
from archive.python_reference_implementation.network_python import (
    compute_msn_python,
    compute_tcs_python,
    compute_mjn_python,
)
```

See main README for current user-facing workflow.
