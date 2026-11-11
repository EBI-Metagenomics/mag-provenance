"""Tests for mag_source_duplicates.py.

ENA is replaced by FakeClient, which answers from the small made-up tables below, so the
tests need no network. Run them from the repo folder with either of:

    python3 -m unittest discover tests
    python3 -m pytest tests
"""

import contextlib
import csv
import gzip
import io
import os
import re
import shutil
import sys
import tempfile
import unittest
import urllib.error
import urllib.parse
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mag_source_duplicates as m  # noqa: E402

# ---------------------------------------------------------------------------
# A fake ENA
#
# FakeClient answers the script's ENA searches from the small tables below, so the tests
# run without a network, in a fraction of a second, and on data that changes in ENA
# can't break them. The accessions are made up for the tests; a few happen to exist in
# ENA, but the tests never contact ENA.
#
# The tables hold one source sample, SAMEA9000001, with the runs ERR9000001 and
# ERR9000002, and samples that point to the source sample in each way the script must
# handle: by its BioSample, its older accession, one of its runs, or through a bin, and
# alongside a dead end, an accession not in ENA, or "missing".
# ---------------------------------------------------------------------------

# A dead end is a sample with no runs and nothing further back.
SAMPLES = {  # BioSample: (older accession, related_sample_accession)
    "SAMEA9000001": ("ERS9000001", ""),                           # the source sample
    "SAMEA9000002": ("ERS9000002", ""),                           # a dead end
    "SAMEA9100001": ("ERS9100001", "SAMEA9000001:derived_from"),  # names the source sample
    "SAMEA9100002": ("ERS9100002", "ERR9000002:derived_from"),    # names one of its runs
    "SAMEA9100003": ("ERS9100003",                                # source + a dead end
                     "SAMEA9000001:derived_from;SAMEA9000002:derived_from"),
    "SAMEA9100004": ("ERS9100004", "SAMEA9000002:derived_from"),  # names only a dead end
    "SAMEA9100005": ("ERS9100005", "ERS9000099:derived_from"),    # names an ERS not in ENA
    "SAMEA9100006": ("ERS9100006", "missing:derived_from"),       # names no accession
    "SAMEA9200001": ("ERS9200001", "ERS9000001:derived_from"),    # names the ERS accession
    "SAMEA9200002": ("ERS9200002", "SAMEA9000001:derived_from"),  # a bin, no genome
    "SAMEA9200003": ("ERS9200003", "SAMEA9200002:derived_from"),  # names the bin
    "SAMEA9200004": ("ERS9200004", "ERR9000001:derived_from"),    # names the other run
    "SAMEA9200005": ("ERS9200005", "SAMEA9000001:derived_from"),  # a non-MAG and a MAG
    "SAMEA9200006": ("ERS9200006", "SAMEA9000001:derived_from"),  # a co-assembly, no MAG
    "SAMEA9200007": ("ERS9200007", "SAMEA9200006:derived_from"),  # names the co-assembly
}
RUNS = [("ERR9000001", "SAMEA9000001", "ERS9000001"),
        ("ERR9000002", "SAMEA9000001", "ERS9000001")]
MAG = "Metagenome-Assembled Genome (MAG)"
# accession: (version, study, sample, assembly_type, run_accession)
ASSEMBLIES = {
    "GCA_900000001": ("1", "PRJEB1", "SAMEA9100001", MAG, "ERR9000001"),
    "GCA_900000002": ("1", "PRJEB2", "SAMEA9100002", MAG, "ERR9000001"),
    "GCA_900000003": ("1", "PRJEB3", "SAMEA9100003", MAG, ""),
    "GCA_900000004": ("1", "PRJEB3", "SAMEA9100004", MAG, ""),
    "GCA_900000005": ("1", "PRJEB4", "SAMEA9100005", MAG, ""),
    "GCA_900000006": ("1", "PRJEB4", "SAMEA9100006", MAG, ""),
    "GCA_900000007": ("29", "PRJNA1", "", "", ""),                 # has no sample
    "GCA_900000008": ("1", "PRJEB4", "SAMEA9100006", "clone or isolate", ""),
    "GCA_900000009": ("1", "PRJEB4", "SAMEA9100006", "COVID-19 outbreak", ""),
    "GCA_900000011": ("1", "PRJEB5", "SAMEA9200001", MAG, "ERR9000002"),
    "GCA_900000012": ("1", "PRJEB6", "SAMEA9200003", MAG, "ERR9000001"),
    "GCA_900000014": ("1", "PRJEB8", "SAMEA9200005", "primary metagenome", ""),
    "GCA_900000015": ("1", "PRJEB8", "SAMEA9200005", MAG, ""),
    "GCA_900000017": ("1", "PRJEB9", "SAMEA9200007", MAG, "ERR9000002"),
}
ANALYSES = {
    "ERZ9000011": ("", "PRJEB5", "SAMEA9200001", MAG, "ERR9000002"),   # GCA_900000011's ERZ
    "ERZ9000013": ("", "PRJEB7", "SAMEA9200004", MAG, "ERR9000001"),
    "ERZ9000016": ("", "PRJEB9", "SAMEA9200006", "primary metagenome", ""),
}
# WGS set: its GCA. Sets with a 4-letter prefix are the older format.
WGS_SETS = {
    "ABCD01000000": "GCA_900000001",     # 4-letter prefix
    "ABCDEF010000000": "GCA_900000003",  # 6-letter prefix
}


def _terms(query):
    """The field="value" pairs of an ENA query:
    'sample_accession="X"' -> [("sample_accession", "X")]."""
    return re.findall(r'(\w+)="([^"]*)"', query or "")


class FakeClient:
    """Answers EnaClient.search() from the tables above, like ENA would."""

    batch_size = 500

    def __init__(self, cache_dir=None):
        self.request_count = 0
        self.cache_hits = 0

    def search(self, result, fields, query=None, accessions=None):
        self.request_count += 1
        wanted = set(accessions or [])
        terms = _terms(query)
        if result == "sample":
            rows = []
            for acc, (sec, rel) in SAMPLES.items():
                links = {x for x in rel.split(";") if x}
                named = any(k == "related_sample_accession" and v in links
                            for k, v in terms)
                hit = acc in wanted or sec in wanted or named
                if hit:
                    rows.append({"sample_accession": acc, "secondary_sample_accession": sec,
                                 "related_sample_accession": rel})
            return rows
        if result == "read_run":
            return [{"run_accession": r, "sample_accession": a,
                     "secondary_sample_accession": s}
                    for r, a, s in RUNS if wanted & {r, a, s}]
        if result in ("assembly", "analysis"):
            table = ASSEMBLIES if result == "assembly" else ANALYSES
            key = "accession" if result == "assembly" else "analysis_accession"
            samples = {v for k, v in terms if k == "sample_accession"}
            rows = []
            for acc, (ver, study, sample, atype, runs) in table.items():
                if acc in wanted or (sample and sample in samples):
                    row = {key: acc, "study_accession": study, "sample_accession": sample,
                           "assembly_type": atype, "run_accession": runs,
                           "center_name": "C_" + study}
                    if result == "assembly":
                        row["version"] = ver
                    rows.append(row)
            return rows
        if result == "wgs_set":
            # Like ENA, reject a 4-letter WGS set in an accession list (HTTP 400 for the
            # whole request), but find it with a query. A 4-letter set is 12 characters,
            # e.g. ABCD01000000.
            if any(len(a) == 12 for a in wanted):
                raise m.EnaError("HTTP 400: accession(s) not valid for wgs_set")
            named = {v for k, v in terms if k == "accession"}
            return [{"accession": acc, "assembly_accession": gca}
                    for acc, gca in WGS_SETS.items() if acc in wanted or acc in named]
        raise AssertionError(result)


class EbiFakeClient(FakeClient):
    """Like FakeClient, but PRJEB4 and PRJEB5 are submitted by EBI (center_name EMG)."""

    def search(self, *args, **kwargs):
        rows = super().search(*args, **kwargs)
        for row in rows:
            if row.get("study_accession") in ("PRJEB4", "PRJEB5"):
                row["center_name"] = "EMG"
        return rows


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_main(args, client=FakeClient):
    """Run the script with a fake ENA and return (exit code, screen output)."""
    screen = io.StringIO()
    with mock.patch.object(m, "EnaClient", client), contextlib.redirect_stderr(screen):
        code = m.main(args)
    return code, screen.getvalue()


def write_lines(path, lines):
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")


class TempDirTest(unittest.TestCase):
    """Gives each test class its own folder, removed at the end."""

    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir)


# ---------------------------------------------------------------------------
# Phase 1 and phase 2, end to end
# ---------------------------------------------------------------------------

class TestBothPhases(TempDirTest):
    """Twelve inputs, covering every source_sample_status and the phase 2 outcomes. The
    clusters file has one id without a version (GCA_900000002)."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        accessions = os.path.join(cls.dir, "acc.txt")
        write_lines(accessions, [
            "GCA_900000001", "GCA_900000002", "GCA_900000003", "GCA_900000004",
            "GCA_900000005", "GCA_900000006", "GCA_900000007", "GCA_900000008",
            "GCA_900000009", "GCA_900000010", "GCF_000001405", "GCA_900000001.1"])
        clusters = os.path.join(cls.dir, "run_clusters.csv")
        write_lines(clusters, [
            "genome_id,cluster", "GCA_900000001.1,c1", "GCA_900000002,c1",
            "GCA_900000011.1,c2", "GCA_900000012.1,c2",
            "GCA_900000015.1,c3", "ERZ9000013,c3"])
        cls.code, cls.screen = run_main(["--accessions", accessions,
                                         "--gemsparcl-clusters", clusters,
                                         "--outdir", cls.dir])
        cls.rows = read_tsv(os.path.join(cls.dir, "mags.tsv"))
        cls.by_input = {r["input_accession"]: r for r in cls.rows if r["input_accession"]}
        cls.by_genome = {r["genome_accession"]: r for r in cls.rows}

    def test_exit_code(self):
        self.assertEqual(self.code, 0)

    def test_source_sample_status_of_each_input(self):
        expected = {
            "GCA_900000001": "traced", "GCA_900000002": "traced",
            "GCA_900000003": "partly_traced",      # names the source sample and a dead end
            "GCA_900000004": "source_not_found",   # names only a dead end
            "GCA_900000005": "source_not_found",   # names a sample not in ENA
            "GCA_900000006": "no_derived_from",    # derived_from is "missing"
            "GCA_900000007": "mag_sample_not_found",
            "GCA_900000008": "not_a_mag", "GCA_900000009": "not_a_mag",
            "GCA_900000010": "accession_not_found",
            "GCF_000001405": "invalid_accession",
            "GCA_900000001.1": "repeated_input",
        }
        for accession, status in expected.items():
            with self.subTest(accession):
                self.assertEqual(self.by_input[accession]["source_sample_status"], status)

    def test_source_group_members(self):
        group = m.fingerprint(["SAMEA9000001"])
        members = sorted(r["genome_accession"] for r in self.rows
                         if r["source_group_id"] == group)
        self.assertEqual(members, sorted([
            "GCA_900000001.1", "GCA_900000002.1", "GCA_900000011.1", "GCA_900000012.1",
            "ERZ9000013", "GCA_900000015.1", "GCA_900000017.1"]))

    def test_erz_copy_of_a_gca_is_dropped(self):
        self.assertNotIn("ERZ9000011", self.by_genome)

    def test_found_non_mag_is_listed(self):
        self.assertEqual(self.by_genome["GCA_900000014.1"]["source_sample_status"],
                         "not_a_mag")

    def test_duplicate_status_after_phase_2(self):
        expected = {
            "GCA_900000001.1": "duplicate", "GCA_900000002.1": "duplicate",
            "GCA_900000011.1": "duplicate", "GCA_900000012.1": "duplicate",
            "GCA_900000015.1": "duplicate", "ERZ9000013": "duplicate",
            "GCA_900000017.1": "no_cluster", "GCA_900000003.1": "no_source_group",
        }
        for genome, status in expected.items():
            with self.subTest(genome):
                self.assertEqual(self.by_genome[genome]["duplicate_status"], status)

    def test_repeated_input_keeps_its_version(self):
        self.assertEqual(self.by_input["GCA_900000001.1"]["genome_accession"],
                         "GCA_900000001.1")

    def test_cluster_matched_without_a_version(self):
        self.assertEqual(self.by_genome["GCA_900000002.1"]["gemsparcl_cluster"], "c1")
        self.assertIn("matched without a version", self.screen)

    def test_no_sequence_check_without_fasta_dir(self):
        self.assertEqual({r["sequence_check"] + r["sequence_md5"] for r in self.rows}, {""})
        self.assertFalse(os.path.exists(os.path.join(self.dir, "identical_sequences.tsv")))

    def test_duplicates_tsv(self):
        dups = read_tsv(os.path.join(self.dir, "duplicates.tsv"))
        runs = {d["gemsparcl_cluster"]: d["declared_runs"] for d in dups}
        self.assertEqual(runs, {"c1": "same", "c2": "different", "c3": "unknown"})
        self.assertIn("center_names", dups[0])

    def test_candidates_after_phase_2(self):
        candidates = read_tsv(os.path.join(self.dir, "organism_check_candidates.tsv"))
        self.assertEqual([c["genome_accession"] for c in candidates], ["GCA_900000017.1"])


# ---------------------------------------------------------------------------
# Reverse lookup from a MAG that names one run
# ---------------------------------------------------------------------------

class TestReverseLookupFromARun(TempDirTest):
    """GCA_900000002's sample names a run (ERR9000002), not the source sample. Other MAGs
    from the same source sample may name the sample, its older accession or its other
    run instead, so the reverse lookup searches with all of them and finds all three."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        accessions = os.path.join(cls.dir, "acc.txt")
        write_lines(accessions, ["GCA_900000002"])
        cls.code, _ = run_main(["--accessions", accessions, "--outdir", cls.dir])
        rows = read_tsv(os.path.join(cls.dir, "mags.tsv"))
        cls.found = {r["genome_accession"] for r in rows if not r["input_accession"]}

    def test_exit_code(self):
        self.assertEqual(self.code, 0)

    def test_finds_the_genome_naming_the_other_run(self):
        self.assertIn("ERZ9000013", self.found)            # names ERR9000001

    def test_finds_the_genome_naming_the_sample(self):
        self.assertIn("GCA_900000001.1", self.found)       # names SAMEA9000001

    def test_finds_the_genome_naming_the_older_accession(self):
        self.assertIn("GCA_900000011.1", self.found)       # names ERS9000001


# ---------------------------------------------------------------------------
# WGS sets as input
# ---------------------------------------------------------------------------

class TestWgsSets(TempDirTest):
    """ENA rejects a WGS set with a 4-letter prefix (e.g. NVUE02000000) in an accession
    list, with HTTP 400 for the whole request, so the script looks WGS sets up with a
    query. Otherwise one such line would stop the whole run. Both prefix lengths must be found and replaced by their GCA, and a WGS set ENA
    doesn't have must get accession_not_found, not stop the run.

        ABCD01     4-letter prefix, short form   -> GCA_900000001
        ABCDEF01   6-letter prefix, short form   -> GCA_900000003
        WXYZ01     not in ENA                    -> accession_not_found
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        accessions = os.path.join(cls.dir, "acc.txt")
        write_lines(accessions, ["ABCD01", "ABCDEF01", "WXYZ01"])
        cls.code, cls.screen = run_main(["--accessions", accessions, "--outdir", cls.dir,
                                         "--no-expand"])
        rows = read_tsv(os.path.join(cls.dir, "mags.tsv"))
        cls.by_input = {r["input_accession"]: r for r in rows}

    def test_exit_code(self):
        self.assertEqual(self.code, 0, self.screen)

    def test_both_prefix_lengths_are_replaced_by_their_gca(self):
        self.assertEqual(self.by_input["ABCD01"]["genome_accession"], "GCA_900000001.1")
        self.assertEqual(self.by_input["ABCDEF01"]["genome_accession"], "GCA_900000003.1")

    def test_wgs_set_not_in_ena(self):
        self.assertEqual(self.by_input["WXYZ01"]["source_sample_status"],
                         "accession_not_found")


# ---------------------------------------------------------------------------
# Identical sequences (--fasta-dir)
# ---------------------------------------------------------------------------

SAME = ">contig_1 MAG\nACGTacgt\nAC\n>contig_2\nTTTT\n"


class TestIdenticalSequences(TempDirTest):
    """Four input genomes have the same sequence. They test the order in which the kept
    copy is chosen: a traced copy first, then one from EMG.

        GCA_900000011                  traced, from EMG       kept
        GCA_900000006                  not traced, from EMG   not kept: traced comes first
        GCA_900000001, GCA_900000002   traced, not from EMG   not kept
    """

    @classmethod
    def write_fasta(cls, name, text):
        path = os.path.join(cls.fasta, name)
        opener = gzip.open if name.endswith(".gz") else open
        with opener(path, "wt") as handle:
            handle.write(text)
        return path

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.fasta = os.path.join(cls.dir, "fasta")
        os.makedirs(cls.fasta)
        # GCA_900000002's file has the same sequence as SAME, written differently:
        # gzipped, with other headers, contig order, line breaks and case. Its checksum
        # must still match. GCA_900000003 is an input with no file in the folder, so it
        # gets no_fasta_file.
        cls.f1 = cls.write_fasta("GCA_900000001.1.fna", SAME)
        cls.f2 = cls.write_fasta("GCA_900000002.fa.gz", ">y\nTTTT\n>x\nACGTA\nCGTAC\n")
        cls.write_fasta("GCA_900000011.1.fasta", SAME)
        cls.write_fasta("GCA_900000006.1.fa", SAME)
        cls.f12 = cls.write_fasta("GCA_900000012.1.fna", ">a\nACGTACGTAA\n>b\nTTTT\n")
        cls.write_fasta("GCA_900000004.2.fna", ">a\nACGT\n")        # another version
        cls.write_fasta("GCA_900000005.fna", "")                      # empty
        cls.write_fasta("GCA_900000007.29.fna", ">a\nGGGGCCCC\n")    # MAG with no sample
        cls.write_fasta("notes.txt", "not a FASTA file")              # ignored

        cls.accessions = os.path.join(cls.dir, "acc.txt")
        write_lines(cls.accessions, [
            "GCA_900000002", "GCA_900000001", "GCA_900000011", "GCA_900000006",
            "GCA_900000012", "GCA_900000004", "GCA_900000005", "GCA_900000003",
            "GCA_900000007"])
        clusters = os.path.join(cls.dir, "run_clusters.csv")
        write_lines(clusters, ["genome_id,cluster", "GCA_900000011.1,c1",
                               "GCA_900000012.1,c1", "GCA_900000001.1,c1",
                               "GCA_900000002.1,c1"])
        cls.out = os.path.join(cls.dir, "res")
        cls.code, cls.screen = run_main(
            ["--accessions", cls.accessions, "--fasta-dir", cls.fasta,
             "--gemsparcl-clusters", clusters, "--outdir", cls.out], client=EbiFakeClient)
        rows = read_tsv(os.path.join(cls.out, "mags.tsv"))
        cls.rows = rows
        cls.by_input = {r["input_accession"]: r for r in rows if r["input_accession"]}

    def column(self, accession, name):
        return self.by_input[accession][name]

    def test_exit_code(self):
        self.assertEqual(self.code, 0)

    def test_checksum_ignores_headers_order_line_breaks_case_and_gzip(self):
        self.assertEqual(m.sequence_checksum(self.f1), m.sequence_checksum(self.f2))

    def test_checksum_differs_for_a_different_sequence(self):
        self.assertNotEqual(m.sequence_checksum(self.f1), m.sequence_checksum(self.f12))

    def test_file_names(self):
        self.assertEqual(m._accession_of_file("GCA_900000001.1.fna.gz"),
                         ("GCA_900000001", "1"))
        self.assertEqual(m._accession_of_file("ERZ9000013.fasta"), ("ERZ9000013", ""))
        self.assertIsNone(m._accession_of_file("notes.txt"))

    def test_sequence_check(self):
        self.assertEqual({a: self.column(a, "sequence_check") for a in self.by_input}, {
            "GCA_900000011": "duplicate_kept", "GCA_900000001": "duplicate",
            "GCA_900000002": "duplicate", "GCA_900000006": "duplicate",
            "GCA_900000012": "unique", "GCA_900000007": "unique",
            "GCA_900000004": "no_fasta_file", "GCA_900000005": "no_fasta_file",
            "GCA_900000003": "no_fasta_file"})

    def test_sequence_md5_for_every_checked_genome(self):
        self.assertEqual(self.column("GCA_900000001", "sequence_md5"),
                         m.sequence_checksum(self.f1))
        self.assertEqual(self.column("GCA_900000012", "sequence_md5"),
                         m.sequence_checksum(self.f12))
        self.assertEqual(self.column("GCA_900000004", "sequence_md5"), "")
        self.assertEqual(self.column("GCA_900000003", "sequence_md5"), "")

    def test_mag_with_no_ena_sample_is_still_checked(self):
        self.assertEqual(self.column("GCA_900000007", "source_sample_status"),
                         "mag_sample_not_found")

    def test_copies_are_still_traced(self):
        self.assertEqual([self.column(a, "source_sample_status") for a in
                          ("GCA_900000001", "GCA_900000002", "GCA_900000006")],
                         ["traced", "traced", "no_derived_from"])
        self.assertEqual(self.column("GCA_900000001", "source_samples"), "SAMEA9000001")

    def test_copies_are_not_compared(self):
        for accession in ("GCA_900000001", "GCA_900000002", "GCA_900000006"):
            with self.subTest(accession):
                self.assertEqual(self.column(accession, "duplicate_status"),
                                 "duplicate_sequence")
                self.assertEqual(self.column(accession, "source_group_id"), "")

    def test_kept_copy_is_compared(self):
        self.assertEqual(self.column("GCA_900000011", "duplicate_status"), "duplicate")

    def test_copies_are_not_added_back_by_reverse_lookup(self):
        self.assertEqual(sum(1 for r in self.rows
                             if r["genome_accession"] == "GCA_900000002.1"), 1)

    def test_copies_are_not_candidates(self):
        candidates = read_tsv(os.path.join(self.out, "organism_check_candidates.tsv"))
        self.assertFalse({"GCA_900000001.1", "GCA_900000002.1"}
                         & {c["genome_accession"] for c in candidates})

    def test_identical_sequences_tsv(self):
        self.assertEqual(read_tsv(os.path.join(self.out, "identical_sequences.tsv")), [{
            "sequence_md5": m.sequence_checksum(self.f1),
            "kept_genome": "GCA_900000011.1", "n_genomes": "4",
            "genomes": "GCA_900000011.1,GCA_900000001.1,GCA_900000002.1,GCA_900000006.1",
            "studies": "PRJEB5,PRJEB1,PRJEB2,PRJEB4"}])

    def test_summary(self):
        for pattern in (r"Step 1: compare the input genomes' sequences",
                        r"\n    unique +2\n", r"duplicate_kept +1  \(one per set",
                        r"\n    duplicate +3\n", r"no_fasta_file +3\n",
                        r"duplicate_sequence +3 +0"):
            with self.subTest(pattern):
                self.assertRegex(self.screen, pattern)

    def test_warnings(self):
        self.assertIn("GCA_900000004.2.fna", self.screen)               # another version
        self.assertIn("GCA_900000005.fna holds no sequence", self.screen)

    def test_missing_folder_stops_the_run(self):
        code, _ = run_main(["--accessions", self.accessions, "--fasta-dir",
                            os.path.join(self.dir, "nope"),
                            "--outdir", os.path.join(self.dir, "res2")])
        self.assertEqual(code, 1)


# ---------------------------------------------------------------------------
# ENA requests: retries, failures and the cache
# ---------------------------------------------------------------------------

class _Answer:
    """A successful HTTP answer, as urlopen() returns it."""

    def __init__(self, text):
        self.text = text

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def read(self):
        return self.text.encode("utf-8")


def http_error(code):
    return urllib.error.HTTPError(m.ENA_SEARCH_URL, code, "error", {},
                                  io.BytesIO(b"ENA's message"))


class ScriptedNetwork:
    """Stands in for urlopen(): each call gets the next outcome from the list, either the
    text of an answer or an exception to raise. Records what was sent."""

    def __init__(self, outcomes):
        self.outcomes = list(outcomes)
        self.sent = []

    def __call__(self, request, timeout=None):
        self.sent.append(urllib.parse.parse_qs(request.data.decode("utf-8")))
        outcome = self.outcomes.pop(0)
        if isinstance(outcome, Exception):
            raise outcome
        return _Answer(outcome)


ROWS = "sample_accession\nSAMEA2622837\n"   # ENA's answer to the test search: one sample


class TestEnaRequests(unittest.TestCase):

    def ask(self, outcomes, client=None, **search):
        """Send one search through a scripted network. Returns the rows, and the network,
        which records what was sent."""
        network = ScriptedNetwork(outcomes)
        client = client or m.EnaClient(retry_wait=0)
        search = search or {"query": 'sample_accession="SAMEA2622837"'}
        with mock.patch.object(m.urllib.request, "urlopen", network):
            rows = client.search("sample", ["sample_accession"], **search)
        return rows, network

    def test_answer_is_read(self):
        rows, network = self.ask([ROWS])
        self.assertEqual(rows, [{"sample_accession": "SAMEA2622837"}])
        self.assertEqual(len(network.sent), 1)

    def test_wrong_request_is_not_retried(self):
        for code in (400, 404, 414):
            with self.subTest(code):
                with self.assertRaisesRegex(m.EnaError, "ENA rejected the request"):
                    self.ask([http_error(code), ROWS])

    def test_server_error_is_retried(self):
        rows, network = self.ask([http_error(503), ROWS])
        self.assertEqual(len(rows), 1)
        self.assertEqual(len(network.sent), 2)

    def test_error_inside_a_normal_answer_is_retried(self):
        rows, network = self.ask(["ERROR occurred while searching", ROWS])
        self.assertEqual(len(rows), 1)
        self.assertEqual(len(network.sent), 2)

    def test_network_failure_is_retried(self):
        rows, network = self.ask([OSError("timed out"), ROWS])
        self.assertEqual(len(rows), 1)

    def test_failure_after_every_attempt_stops_the_run(self):
        # So a failure is never read as "ENA found nothing".
        client = m.EnaClient(retry_wait=0)
        for outcome in (http_error(503), "ERROR occurred", OSError("timed out")):
            with self.subTest(repr(outcome)):
                network = ScriptedNetwork([outcome] * client.max_retries)
                with mock.patch.object(m.urllib.request, "urlopen", network):
                    with self.assertRaisesRegex(m.EnaError, "failed after 4 attempts"):
                        client.search("sample", ["sample_accession"], query="x")
                self.assertEqual(len(network.sent), client.max_retries)

    def test_every_attempt_is_counted(self):
        client = m.EnaClient(retry_wait=0)
        self.ask([http_error(500), ROWS], client=client)
        self.assertEqual(client.request_count, 2)

    def test_cache_answers_the_same_request_again(self):
        with tempfile.TemporaryDirectory() as cache:
            client = m.EnaClient(cache_dir=cache, retry_wait=0)
            self.ask([ROWS], client=client)
            # The network has no answer left, so this one must come from the cache.
            rows, network = self.ask([], client=client)
            self.assertEqual(rows, [{"sample_accession": "SAMEA2622837"}])
            self.assertEqual(client.cache_hits, 1)
            self.assertEqual([f for f in os.listdir(cache) if f.endswith(".part")], [])

    def test_failed_answer_is_not_cached(self):
        with tempfile.TemporaryDirectory() as cache:
            client = m.EnaClient(cache_dir=cache, retry_wait=0, max_retries=1)
            with self.assertRaises(m.EnaError):
                self.ask(["ERROR occurred"], client=client)
            self.assertEqual(os.listdir(cache), [])

    def test_sample_accession_forms_go_in_separate_requests(self):
        _, network = self.ask(["", ""], accessions=["SAMEA2622837", "ERS811241"])
        sent = sorted(s["includeAccessions"][0] for s in network.sent)
        self.assertEqual(sent, ["ERS811241", "SAMEA2622837"])

    def test_accessions_are_sent_in_batches(self):
        client = m.EnaClient(retry_wait=0, batch_size=2)
        accessions = ["GCA_90000000%s" % i for i in range(1, 6)]
        _, network = self.ask(["", "", ""], client=client, accessions=accessions)
        self.assertEqual([s["includeAccessions"][0].count(",") + 1 for s in network.sent],
                         [2, 2, 1])

    def test_ena_failure_ends_the_run_with_an_error(self):
        class FailingClient(FakeClient):
            def search(self, *args, **kwargs):
                raise m.EnaError("ENA request failed after 4 attempts")

        with tempfile.TemporaryDirectory() as out:
            accessions = os.path.join(out, "acc.txt")
            write_lines(accessions, ["GCA_900000001"])
            code, screen = run_main(["--accessions", accessions, "--outdir", out],
                                    client=FailingClient)
            self.assertEqual(code, 1)
            self.assertIn("ERROR: ENA request failed", screen)
            self.assertFalse(os.path.exists(os.path.join(out, "mags.tsv")))


if __name__ == "__main__":
    unittest.main()
