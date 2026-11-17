# `mags.tsv`: columns and statuses

`mags.tsv` has one row per genome: every accession in your list (an exact repeat of a line
is skipped), and every genome found by the reverse lookup. Phase 1 writes it, and phase 2
updates it. The [README](../README.md) explains how the script works.

## Columns

Copied from ENA: the genome's accession (`genome_accession`), its sample (`mag_sample`),
the runs it lists (`declared_runs`), study, center, last update, assembly and binning
software, the submitter's own completeness and contamination scores, assembly type,
taxonomy and FASTA link (`fasta_ftp`). `genome_accession` is written with its version
(GCA_937896655.1), which phase 2 uses to match the clusters file.

Added by the script:

| Column | |
|---|---|
| `input_accession` | the accession as written in your accessions file; empty for genomes found by the reverse lookup |
| `source_samples` | the source samples the genome traces back to |
| `source_group_id` | its source group, e.g. `sg_3bd7455b797d`; the same samples give the same id in every run |
| `gemsparcl_cluster` | its cluster, after phase 2; empty if the clusters file doesn't have it |
| `source_sample_status` | whether it could be traced (below) |
| `duplicate_status` | what the comparison found (below) |
| `sequence_check` | the result of the sequence comparison (below) |
| `sequence_md5` | the checksum of the genome's sequences, for matching against later runs |

## `source_sample_status`

| Value | |
|---|---|
| **`traced`** | every `derived_from` link led to a source sample |
| `partly_traced` | some links led to source samples, but at least one didn't, so the genome gets no source group |
| `mag_sample_not_found` | ENA lists no sample for the genome, or has no record of that sample (for example, reference genomes such as GCA_000001405, or NCBI samples not yet in ENA, such as SAMN30429671) |
| `no_derived_from` | the genome's sample has an empty `derived_from`, or one with no accession (for example, "missing") |
| `source_not_found` | `derived_from` names an accession, but no sample with runs was found behind it: ENA has no record of that accession (for example, ERS950426), or it is a sample with no runs and nothing further back (for example, SRS1615213) |
| `trace_too_long` | no sample with runs within 10 steps |
| `not_a_mag` | ENA's assembly type says it isn't a MAG (see [Limits](../README.md#limits)) |
| `repeated_input` | the same genome as an earlier input line, for example its WGS set or another version |
| `accession_not_found` | a GCA, ERZ or WGS accession ENA has no record of |
| `invalid_accession` | not a GCA, ERZ or WGS accession |

## `duplicate_status`

Phase 1 sets `to_compare`, `single_study`, `duplicate_sequence` or `no_source_group`.
Phase 2 replaces every `to_compare` with one of its four values; the other phase 1 values
stay.

| Value | Set by | |
|---|---|---|
| **`duplicate`** | phase 2 | a genome from another study in its source group has the same gemsparcl cluster |
| `not_duplicate` | phase 2 | the genomes from other studies in its source group all have other clusters |
| `other_study_no_cluster` | phase 2 | no match, but a genome from another study has no cluster, so it could still be a duplicate |
| `no_cluster` | phase 2 | its source group is used by 2 or more studies, but this genome isn't in the clusters file |
| `to_compare` | phase 1 | its source group is used by 2 or more studies, so phase 2 compares it |
| `single_study` | phase 1 | all genomes in its source group are from one study, so it can't be a duplicate |
| `duplicate_sequence` | phase 1 | it has the same sequence as another input genome, which is compared instead |
| `no_source_group` | phase 1 | its `source_sample_status` isn't `traced`, so it can't be compared |

## `sequence_check`

| Value | |
|---|---|
| `unique` | no other input genome has the same sequence |
| `duplicate_kept` | one of a set of input genomes with the same sequence: the copy kept for the comparisons |
| `duplicate` | another copy of that set; it gets `duplicate_status` `duplicate_sequence` |
| `no_fasta_file` | no usable FASTA file: missing, another version, empty or unreadable |
| (empty) | not checked: no `--fasta-dir`, a genome found by the reverse lookup, or an input that is `not_a_mag`, `repeated_input`, `accession_not_found` or `invalid_accession` |

The kept copy is chosen in this order: a traced genome first (only traced genomes get a
source group), then one submitted by EMG, then a GCA before an ERZ, then the lowest
accession.
