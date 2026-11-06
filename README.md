# mag-provenance

Finds genomes in ENA that are copies of one another: the same organism, built from the BioSample, but assembled and deposited separately in different projects.

A MAG (metagenome-assembled genome) is a genome reconstructed from a mixed environmental
sample. The raw sequencing data behind it stays public in ENA,
so other groups often download the same reads, run their own assembly and binning, and
deposit the genomes they get.

The archive then holds several versions of the same organism from the same original sample, counting them separately overstates how
much has really been sampled.

---

## What this tool treats as a duplicate

| | | How we know |
|---|---|---|
| 1 | The genomes come from the same original sample | same BioSample|
| 2 | Different groups made them | they were deposited under different study accessions |
| 3 | They are the same organism | same gemsparcl cluster |

This tool works out 1 and 2 by reading ENA, in two steps.

**Walking up.** A MAG is not filed under the sample it came from. It gets a sample record
of its own, which names what it was derived from in a field called `derived_from`. The tool follows that field, sometimes through more than one step, until it arrives at a **BioSample** with sequencing runs attached. That is a source of the genome. `derived_from` can name several samples at once, because a genome may have been co-assembled
from more than one of them together. Each one is followed separately, and they can be
different distances away, so a genome can end up with several source BioSamples. Each MAG is grouped by an identical source-sample set. The walk-up ends to an empty `derived_from`. If nothing at all was found, the genome is reported as untraceable.

**Reverse look-up.** It then asks ENA which *other* genomes were built from those same
samples. 

For 3 it uses a clustering file, made by
gemsparcl.
This script never runs gemsparcl and never downloads a genome.

---

## Requirements

Python 3.7 or newer. The script needs
to reach the ENA website while running.

---

## Running it

```bash
python3 mag_source_duplicates.py \
    --accessions accessions.txt \
    --gemsparcl-results /path/to/gemsparcl_output/ \
    --outdir results \
    --cache-dir cache
```

| Option | |
|---|---|
| `--accessions` | your list of genome accessions, one per line |
| `--gemsparcl-results` | folder holding gemsparcl's output files (`*_clusters.csv`, `*.dists`) |
| `--outdir` | where to write the results (default: `results`) |
| `--cache-dir` | keeps the ENA replies on disk, so a run that is interrupted can restart cheaply instead of asking again |
| `--no-expand` | only look at the genomes you listed. Much faster, but finds far fewer duplicates — see below |

Leaving out `--gemsparcl-results` stops after the ENA part. You still get the groups of
genomes that share a sample, and a list of which ones are worth comparing, but no genome is
called a duplicate because nothing has checked whether they are the same organism.

### The input file

One accession per line. Blank lines, `#` comments and extra columns are ignored, so an NCBI
`assembly_summary` file can be used as it is.

```
GCA_937896655
ERZ879826
CABKLP010000000      # a WGS set accession also works
```

Three kinds of accession are accepted: `GCA_` assemblies, `ERZ` analyses, and WGS sets.
Anything else is reported in the output as `invalid_accession`; it does not stop the run.

---

## What you get

| File | |
|---|---|
| **`duplicates.tsv`** | one row per set of genomes that are copies of each other |
| `mags.tsv` | one row per genome, including the ones that could not be traced |
| `source_groups.tsv` | one row per set of genomes that came from the same sample |
| `organism_check_candidates.tsv` | genomes not present in the gemsparcl clustering that need to be queried separately by gemsparcl, with a link to download each |

In `mags.tsv`, every genome is labelled:

| Label | Meaning |
|---|---|
| **`duplicate_independent`** | same organism, same sample, different studies |
| `duplicate_same_study` | same organism and sample, but one study — redundant work within one group |
| `co_derived_different_organism` | came from the same sample but is a different organism. |
| `organism_unknown` | we could not check the organism.|
| `pending_organism_check` | still waiting to be compared with `gemsparcl query` |
| `same_study_redundancy`, `singleton`, `partial_source`, `unresolved` | see the `resolution_status` column for why |

Every run prints two numbers: how much of your list could be traced back to a sample, and
how many of the samples were used by more than one group.

