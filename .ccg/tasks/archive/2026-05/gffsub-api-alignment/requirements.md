# gffsub API alignment

## User Request

Check whether the refactored embedded `deps/gffsub` API in this repository is
aligned with the reference project at `F:/codex/gffsub`.

## Acceptance Criteria

- Identify public headers/classes/functions exposed by the reference `gffsub`.
- Compare embedded `deps/gffsub` headers and source against the reference.
- Check haplokit C++ call sites that depend on gffsub APIs.
- Report API mismatches, compatibility risks, and concrete next steps.

## Constraints

- Start read-only unless a mismatch is clear and safely fixable.
- Do not touch unrelated pending packaging/release edits.
