# mag-provenance

`mag_source_duplicates.py` finds duplicate MAGs in ENA: genomes from the same source
data, deposited by different studies, that gemsparcl assigns to the same cluster.

The reads behind a MAG stay public in ENA, so another group can download the same reads,
assemble and bin them again, and deposit the genomes under its own study. ENA then holds
several genomes of the same organism from the same sample, and counting each copy
separately overstates how many different genomes were recovered.

---

## Quick start

You need Python 3.7 or newer (no extra packages). Phase 1 needs internet access to the
[ENA Portal API](https://www.ebi.ac.uk/ena/portal/api/). Phase 2 needs a gemsparcl
clusters file; the script doesn't run gemsparcl.

```bash
# phase 1: look the genomes up in ENA; with --fasta-dir (optional), also find genomes
# with the same sequence by comparing checksums of their FASTA files
python3 mag_source_duplicates.py --accessions accessions.txt --fasta-dir fasta/ \
    --outdir results --cache-dir cache

# phase 2: compare the genomes' clusters in a gemsparcl clusters file
python3 mag_source_duplicates.py --gemsparcl-clusters gemsparcl_out/run_clusters.csv \
    --outdir results

# or both phases in one go
python3 mag_source_duplicates.py --accessions accessions.txt --fasta-dir fasta/ \
    --gemsparcl-clusters gemsparcl_out/run_clusters.csv --outdir results --cache-dir cache
```

The main result is `results/duplicates.tsv`, one row per group of duplicates.

---

## What counts as a duplicate

Genomes are duplicates when all three of these hold:

| | | How the script knows |
|---|---|---|
| 1 | They come from the same source data | they trace back to the same source samples in ENA |
| 2 | Different studies deposited them | they have different study accessions |
| 3 | They are the same organism | gemsparcl assigned them to the same cluster |

The **source sample** is the sample that holds the sequencing runs a MAG was built from.
A co-assembled MAG is built from several samples, so it has several source samples. Two
genomes count as built from the same source data only when they have exactly the same
set of source samples.

The script compares samples, not runs: two genomes from SAMEA2620995 count as built from
the same source data even if they used different runs of it (ERR598947, ERR599131 or
ERR10235166). The `declared_runs_match` column of `duplicates.tsv` shows whether the
genomes in a duplicate group list the same runs in ENA.

Genomes deposited by the same study are not counted as duplicates. The exception is
genomes with the same sequence: with `--fasta-dir`, one copy is kept and the others get
`duplicate_sequence`, whatever their study (see `sequence_check` in
[docs/output.md](docs/output.md)).

---

## How it works

The script runs in two phases. Phase 1 searches ENA; phase 2 reads `mags.tsv` from phase 1
and a gemsparcl clusters file that you provide. Run the phases one after the other, or give both options to run
them in one go. The script never runs gemsparcl and never downloads a genome.

### Phase 1: ENA (checks 1 and 2)

**Step 1: compare sequences (with `--fasta-dir`).** The same genome is sometimes
deposited more than once. The script computes a checksum of each input MAG's sequences
from its FASTA file. Genomes with the same checksum are compared only once in phase 2:
one copy is kept, and the others get `duplicate_sequence`. Every copy is still traced in
step 2, so its provenance is recorded.

**Step 2: trace each genome back to its source samples.** A MAG has its own sample in
ENA, and that sample's `derived_from` field names what the MAG was built from. The script
follows `derived_from` until it reaches a sample with runs: that is a source sample. For
example, GCA_937896655 has the sample SAMEA14075323, whose `derived_from` names
SAMEA2620995, a sample with runs. So the source sample of GCA_937896655 is SAMEA2620995.

- The trace can take several steps, such as MAG sample → bin sample → source sample. The
  trace stops after 10 steps (`trace_too_long`).
- When `derived_from` names a run, the run's sample is the source sample.
- When `derived_from` names several samples (a co-assembly), each one is followed, so a
  genome can have several source samples.

Genomes with exactly the same source samples form a **source group**.

**Step 2, continued: the reverse lookup.** Genomes that other studies built from the same
source samples are usually not in your list. So the script searches ENA for every sample
whose `derived_from` names a source sample, by any of its accessions: SAMEA2620995, its
older accession ERS490230, or one of its runs. A sample found this way that has no MAG,
such as a bin, is searched from in the same way, up to 10 steps. Every genome found is
traced like the input genomes, and joins a source group only if it has exactly the same
source samples. Genomes found that aren't MAGs, by their ENA assembly type, are listed as
`not_a_mag` and not traced.

At the end of phase 1, every genome in a source group used by 2 or more studies has
`duplicate_status` `to_compare`.

### Phase 2: gemsparcl clusters (check 3)

Phase 2 reads a gemsparcl clusters file (`*_clusters.csv`, with columns `genome_id` and
`cluster`). Within each source group used by 2 or more studies, genomes from different
studies that are in the same gemsparcl cluster are duplicates. Genomes in those groups
that the clusters file doesn't have are listed in `organism_check_candidates.tsv`, with a
download link, so you can run them through `gemsparcl query`.

---

## Options

| Option | |
|---|---|
| `--accessions FILE` | run phase 1 on this list of genome accessions |
| `--fasta-dir DIR` | phase 1: compare the input MAGs' sequences, using the FASTA files in this folder (see below) |
| `--gemsparcl-clusters FILE` | run phase 2 with this gemsparcl clusters file (columns `genome_id` and `cluster`) |
| `--outdir DIR` | folder for the output files; phase 2 reads `mags.tsv` from here (default: `results`) |
| `--cache-dir DIR` | save the metadata ENA returns for each search in this folder, so a rerun or an interrupted run reuses it instead of searching ENA again |
| `--no-expand` | skip the reverse lookup: only the genomes in your list are compared, so far fewer duplicates are found |

---

## Input files

### The accessions file (for `--accessions`)

A text file with one accession per line. A version is optional: the script always uses
the genome's current version in ENA. Blank lines, `#` comments and extra columns (after a
tab or comma) are ignored.

```
# my MAGs
GCA_937896655.1
ERZ21930820
CALYNH010000000      # a WGS set; replaced by its GCA, GCA_945214255
```

Three kinds of accession are accepted: `GCA_` assemblies, `ERZ` analyses and WGS sets.
Anything else, such as `GCF_000001405`, gets `invalid_accession` and doesn't stop the run.
An NCBI `assembly_summary` file can be used as it is: its header lines start with `#`, and
only its first column, the accession, is read.

### The FASTA folder (for `--fasta-dir`)

A folder with one FASTA file per genome, plain or gzipped, named by the genome's
accession:

```
fasta/
    GCA_937896655.1.fna.gz
    GCA_945214255.fasta.gz
    ERZ21930820.fa
```

Files are matched to genomes by name: the accession, an optional version, then `.fa`,
`.fasta`, `.fna` or `.fsa`, optionally followed by `.gz`. For GCA_937896655, whose current
version in ENA is 1:

| File | Used? |
|---|---|
| `GCA_937896655.1.fna.gz` | yes |
| `GCA_937896655.fasta.gz` | yes: a name without a version matches whatever version is current in ENA, so the folder should be recent |
| `GCA_937896655.2.fna.gz` | no: another version; a warning names it |

A WGS set input is matched by its GCA. Other files in the folder are ignored. The
checksum doesn't depend on FASTA headers, contig order, line length or upper and lower
case.

### The gemsparcl clusters file (for `--gemsparcl-clusters`)

A comma-separated file with a header row. Phase 2 reads two columns, `genome_id` (the
genome accession, with its version) and `cluster`; any other columns are ignored. For
example:

```
genome_id,cluster
GCA_937896655.1,683
GCA_945214255.1,1728
ERZ21930820,2537
```

gemsparcl groups genomes of the same organism into the same cluster. A genome that isn't
in the file gets no cluster. Phase 2 also reads `mags.tsv` from `--outdir`, written by phase 1.

---

## Output files

All files are written to `--outdir`.

| File | Written by | |
|---|---|---|
| **`duplicates.tsv`** | phase 2 | one row per group of duplicates |
| `mags.tsv` | phase 1, updated by phase 2 | one row per genome: every accession in your list (an exact repeat of a line is skipped), and every genome found by the reverse lookup |
| `identical_sequences.tsv` | phase 1, with `--fasta-dir` | one row per set of input genomes with the same sequence |
| `organism_check_candidates.tsv` | phase 2 | genomes to run through `gemsparcl query` |
| `mag_source_duplicates.log` | both | what was printed on screen, plus step-by-step details; each run is added to the end |

### `duplicates.tsv`

| Column | |
|---|---|
| `duplicate_id` | an id for the duplicate group (the genomes in this row); the same genomes get the same id in every run |
| `source_group_id` | the source group the genomes share |
| `gemsparcl_cluster` | the cluster the genomes share |
| `n_genomes`, `n_studies` | how many genomes and studies are in the group |
| `declared_runs_match` | whether the runs the genomes list in ENA are the same (below); for information only |
| `genomes`, `studies`, `center_names` | the genomes, their studies and the centers that submitted them |

`declared_runs_match` is `same` when every genome lists runs and the lists are the same,
`different` when they all list runs but not the same ones, and `unknown` when at least
one genome lists none. It doesn't decide anything, because different runs of the same
source sample count as the same source data.

### `mags.tsv`

One row per genome, with its ENA metadata, source samples, source group, gemsparcl
cluster and three status columns: `source_sample_status` (could it be traced?),
`duplicate_status` (what the comparison found) and `sequence_check` (the sequence
comparison). Every column and status value is described in
[docs/output.md](docs/output.md).

### `identical_sequences.tsv`

| Column | |
|---|---|
| `sequence_md5` | the checksum the genomes share |
| `kept_genome` | the copy kept for the comparisons |
| `n_genomes`, `genomes` | how many genomes have this sequence, and which, kept copy first |
| `studies` | their studies |

### `organism_check_candidates.tsv`

Genomes in source groups used by 2 or more studies that have no gemsparcl cluster, so you
can assign them to the existing clusters with `gemsparcl query` after phase 2. Columns:
`source_group_id`, `genome_accession` and `fasta_url` (a download link).

---

## The screen summary

Each phase prints a summary, which is also saved in `mag_source_duplicates.log`. Phase 1
gives the `sequence_check` counts (step 1), the `source_sample_status` counts and the
reverse lookup (step 2), then the `duplicate_status` counts. Phase 2 gives the
`duplicate_status` counts of the genomes it compared. The `duplicate_status` tables count
the genomes in your list and the genomes found by the reverse lookup in separate columns.
For example, the phase 2 counts of a run on 1,000 genomes:

```
  duplicate_status                       in your list   reverse lookup
    duplicate                                      63            2,057
    not_duplicate                                   0                1
    other_study_no_cluster                        101            2,041
    no_cluster                                     78            5,597
    total                                         242            9,696
```

The traced percentage leaves out inputs that can't be traced at all (`not_a_mag`,
`repeated_input`, `accession_not_found`, `invalid_accession`), and names how many were
left out.

---

## Limits

- **What counts as a MAG.** ENA's `assembly_type` has eight values, or is empty:
  - MAGs: `Metagenome-Assembled Genome (MAG)` and `binned metagenome`;
  - not MAGs, so `not_a_mag` and not traced: `primary metagenome`, `clone or isolate`,
    `environmental single-cell amplified genome (SAG)`, `metatranscriptome`,
    `COVID-19 outbreak` and `clinical isolate assembly`.

  For now, genomes with an empty `assembly_type` are also counted as MAGs. Some isolates
  have it empty too (for example, GCA_000008125, Thermus thermophilus HB27), so they are
  traced as well. Whether to keep counting them as MAGs is to be confirmed, or replaced by
  a better way to tell MAGs apart in ENA.
- **Only input MAGs are checked for identical sequences.** The reverse lookup can find far
  more genomes than the input list, and hashing all their files would make every run much
  longer, so the genomes it finds are not hashed. Nor are `not_a_mag` inputs, as this version is
  for MAGs. An identical copy found by the reverse lookup is therefore compared in phase 2
  like any other genome, and can end up in the `gemsparcl query` list.
- **A sample range** in `derived_from` (for example, ERS9142775-ERS9142786) is not
  expanded. Only genomes that list exactly the same range share a source group.
- **A partly traced genome gets no source group.** At least one of its `derived_from`
  links didn't lead to a sample with runs: for example, the link names an accession ENA
  has no record of, or is a value such as "missing". Its full set of source samples is therefore
  unknown, so it can't be compared. The source samples that were found are still used in
  the reverse lookup.
- **A MAG that has both an ERZ analysis and a GCA** under the same sample, such as
  ERZ21930820 and GCA_963687975, is reported once by the reverse lookup, as the GCA.
- **Accessions without a version in the clusters file.** Phase 2 matches each genome
  accession in the clusters file (column `genome_id`) with its version: `GCA_937896655.1`
  matches only `GCA_937896655.1`. If the clusters file has accessions without a version,
  `GCA_937896655` matches the version in `mags.tsv`, and a warning gives the number of
  genomes matched this way.
- **Run phase 1 again after updating the script.** Phase 2 reads the `mags.tsv` that
  phase 1 wrote. If that file comes from an older version of the script with other
  columns, phase 2 stops with an error. If the columns are the same but a status was
  renamed, phase 2 doesn't notice, and the old name stays in `mags.tsv`.

---

## Tests

```bash
python3 -m unittest discover tests
```

The tests fake ENA and use temporary files, so they need no network access and no FASTA
folder.
