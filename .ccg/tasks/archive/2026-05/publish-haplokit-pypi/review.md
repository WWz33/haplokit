# Review

## Result

- PyPI already had `haplokit 0.1.1`, so the release version was bumped to `0.1.2`.
- Built `haplokit-0.1.2.tar.gz` and a local Linux wheel for validation.
- Selected sdist-only upload for PyPI because this package is Linux-source-build oriented and the local wheel is CPython/platform-specific.

## Validation

- `twine check /mnt/f/tmp/haplokit-dist-0.1.2/*` passed.
- Clean Linux sdist install passed.
- `haplokit view data/var.sorted.vcf.gz -r scaffold_1:4300-5000` produced `hap_summary.tsv`.
