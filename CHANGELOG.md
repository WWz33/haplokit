# Changelog

All notable changes to haplokit are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions adhere to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- Missing VCF/BCF index no longer aborts with SIGSEGV. htslib increments the
  reader count before the index check fails, so the error path destroyed a
  half-built reader whose `aux->closefile` array was still NULL; the reader
  count is now reset before destruction and the command exits with the
  `failed to add indexed VCF/BCF reader` message
- TCS intermediate vertices no longer read freed memory. `insert_intermediates`
  cached `const Vertex*` pointers across `add_median_vertex` calls, and the
  `vertices_` vector reallocation invalidated them; the endpoint sequences are
  copied before the loop
- Network distances on pipe-delimited state strings (`0/0|1/1|...`) count one
  difference per variant site. Token-wise comparison replaces character-wise
  comparison, so a single SNP is weight 1 instead of 2, and the SIMD path is
  restricted to pure nucleotide input where 2-bit packing is meaningful
  (previously a 33-site-or-longer state string packed to all-zero words and
  reported distance 0)
- BED left coordinates are treated as 0-based half-open, matching the BED
  specification and the GFF annotation path in the same function; a BED row
  `chr1 4300 4345` now selects 1-based 4301-4345 instead of 4300-4345
- `--max-diff` keeps per-site states in the summary rows. The cluster label is
  no longer emitted as the haplotype pattern, so `hap_summary.tsv` rows match
  the header width, and networks built from `--max-diff` runs compare allele
  patterns instead of cluster identifiers
- `--max-diff 0` is treated as exact grouping. It previously took the
  approximate path, reported `approx-region`, and produced different haplotype
  keys than the same command without the flag
- Haplotype frequency denominators cover the analyzed samples only. Population
  totals came from the whole population file and the overall total from every
  sample in the region, so `--samples-file` reported `wild 1/12` under a header
  of `Individuals: 3`, and samples dropped by strict grouping left frequencies
  that did not sum to 1
- TCS degree-2 intermediate vertices are deleted after their paths are
  collapsed. The collapse removed the edges but kept the vertices, leaving
  every collapsed median in the output as an isolated node (61 of 69 nodes on
  the bundled `scaffold_1:4300-5209` example)
- MJN obsolete median vertices are deleted instead of being left isolated, and
  the iteration terminates when no vertex is removed rather than running to its
  iteration cap
- The network backend exits non-zero on malformed JSON input. It previously
  printed an empty network and returned success, which the Python caller could
  not distinguish from a valid empty result
- Inverted and zero-based regions are rejected with an error message.
  `scaffold_1:5209-4300` returned an empty result and exit code 0 on
  `view-json`, while `view-targets-json` rejected the same input; non-numeric
  coordinates surfaced the raw `stoi` message

### Added

- `--network-epsilon` for `haplokit view`, the tolerance (in mutation steps
  above a triplet's Steiner bound) for adding median vertices to an MJN
  network. Without it the CLI always passed epsilon 0, where no median can
  qualify and `-N mjn` produced the same graph as `-N msn`
- `analyzed_sample_count` in the JSONL payload: the number of samples that
  received a haplotype, i.e. the denominator behind the reported frequencies
- `hap_summary.tsv` `freq` column note: it holds the sample count, matching the
  `Accession` column beside it

### Changed

- Table rendering recognises haplotype rows by excluding the `POS` and
  `ALLELE` metadata rows instead of testing for a leading `H`. A custom
  `--hap-prefix` (for example `-H A`) previously left the sample list in the
  count column, dropped the `n/N` values, and skipped indel abbreviation
- Dead code removed: the duplicate `src/cpp/third_party/json.hpp` (byte-identical
  to `deps/nlohmann/nlohmann/json.hpp`), the uncalled graph helpers
  (`compute_distance_matrix`, `floyd_warshall_subset`, `shortest_path_length`,
  `Graph::clear`), and the unused TCS helpers (`compute_median_base`,
  `update_shortest_paths`), plus `MJNAlgorithm::is_feasible_median`, which was
  declared but never defined
- The no-op `bcf_sr_set_samples` calls were removed. The synced reader only
  subsets records in streaming mode, and `BCF_SR_REQUIRE_IDX` disables
  streaming, so the calls never affected the records returned
