# Requirements

- Fix runtime failures such as `libbz2.so.1.0: cannot open shared object file` after pip source installs in conda/mamba environments.
- Preserve the existing backend discovery and build behavior.
- Apply the runtime environment consistently to `haplokit_cpp` and `haplokit_network_backend` subprocess calls.
- Document the immediate workaround and release the fix as the next PyPI version.
