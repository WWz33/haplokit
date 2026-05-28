# Requirements

Rework CLI help organization so it is easier to scan without renaming haplokit's domain-specific commands.

Scope:
- Keep existing command names and argument compatibility.
- Do not blindly copy external command names.
- Organize `view --help` into haplotype-oriented groups.
- Organize `phenotype stat/box --help` into haplotype/phenotype-oriented groups.
- Preserve the newly added short options.
- Add tests proving grouped help headings remain present.

