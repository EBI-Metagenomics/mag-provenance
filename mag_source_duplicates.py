#!/usr/bin/env python3
"""
Find MAGs in ENA that are the same organism, built from the same source sample, and
deposited more than once.

Public reads can be assembled and binned again by another study, and the genomes
deposited under that study. Two genomes are duplicates when all three of these
hold:

    1. they trace back to the same source sample(s)   (ENA metadata)
    2. they were deposited under different studies    (ENA metadata)
    3. they are the same organism                     (same gemsparcl cluster)

A source sample is the sample that holds the sequencing runs a MAG was built from. A
sample can have several runs (e.g. extra sequencing lanes), and ENA often does not say
which of them a MAG used. So genomes from the same source sample count as coming from
the same source data, even if they used different runs.

The script works in two phases. You can run them one at a time, or both in one go by
giving both options:

    Phase 1 (--accessions) checks 1 and 2. It looks each genome up in ENA, follows its
    sample's `derived_from` links back to the source sample, then asks ENA which other
    genomes came from that sample (the reverse lookup). Genomes with exactly the same
    source samples form a "source group". Writes mags.tsv. With --fasta-dir, input MAGs
    with the same sequence are compared only once (see identical_sequences.tsv).

    Phase 2 (--gemsparcl-clusters) checks 3. Within each source group, genomes that
    gemsparcl put in the same cluster become a duplicate group. Writes duplicates.tsv
    and updates mags.tsv. Genomes in source groups with 2 or more studies that are not
    in the clusters file are listed in organism_check_candidates.tsv, so they can be run
    through `gemsparcl query`.
"""

import argparse
import csv
import gzip
import hashlib
import http.client
import json
import logging
import os
import re
import shlex
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from dataclasses import dataclass, field, fields as dataclass_fields, replace
from typing import List, Set

LOG = logging.getLogger("mag_source_duplicates")
LOG_FILE = "mag_source_duplicates.log"


def _dedupe(items):
    """Remove duplicates but keep the original order, so output is the same on every run."""
    seen = set()
    out = []
    for item in items:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


# ---------------------------------------------------------------------------
# Input accessions (phase 1)
# ---------------------------------------------------------------------------

KIND_GCA = "GCA"
KIND_ERZ = "ERZ"
KIND_WGS = "WGS"
KIND_INVALID = "invalid"

# A WGS set is a 4- or 6-letter prefix, a 2-digit version, then zeros. ENA rejects the
# short form (CALBOS01), so we pad it.
# GCF_ (RefSeq) is not accepted: ENA does not hold RefSeq, and one GCF_ in a batch makes
# ENA reject the whole batch.
RE_GCA = re.compile(r"^GCA_[0-9]{9}$")
RE_ERZ = re.compile(r"^ERZ[0-9]+$")
RE_WGS_PADDED = re.compile(r"^[A-Z]{4}([A-Z]{2})?[0-9]{2}0{6,9}$")
RE_WGS_SHORT = re.compile(r"^[A-Z]{4}([A-Z]{2})?[0-9]{2}$")

# Zeros to add to a short WGS accession, by its length: AAAA02 -> AAAA02000000,
# CALBOS01 -> CALBOS010000000.
WGS_PAD = {6: "0" * 6, 8: "0" * 7}


@dataclass(frozen=True)
class Accession:
    """One input line: what the user wrote (raw), its kind (GCA, ERZ, WGS or invalid), and
    the form sent to ENA (query)."""

    raw: str
    kind: str
    query: str


def classify_accession(raw):
    """Check one line of the accessions file and return the accession to send to ENA.

    Accepted: a GCA assembly, an ERZ analysis, or a WGS set. Examples:

        GCA_937896655.1  -> GCA_937896655    version removed; ENA's latest version is used
        ERZ21930252      -> ERZ21930252
        CALBOS01         -> CALBOS010000000  padded with zeros; ENA rejects the short form

    Anything else, such as a typo or a GCF_ accession, gets kind "invalid". The run carries
    on, and the line appears in mags.tsv as invalid_accession.
    """
    original = raw.strip() if raw else ""
    token = original.upper()
    if not token:
        return Accession(raw=original, kind=KIND_INVALID, query="")

    unversioned = split_version(token)[0]

    if RE_GCA.match(unversioned):
        return Accession(raw=original, kind=KIND_GCA, query=unversioned)
    if RE_ERZ.match(unversioned):
        return Accession(raw=original, kind=KIND_ERZ, query=unversioned)
    if RE_WGS_PADDED.match(unversioned):
        return Accession(raw=original, kind=KIND_WGS, query=unversioned)
    if RE_WGS_SHORT.match(unversioned):
        return Accession(raw=original, kind=KIND_WGS,
                         query=unversioned + WGS_PAD[len(unversioned)])
    return Accession(raw=original, kind=KIND_INVALID, query="")


# ---------------------------------------------------------------------------
# Reading derived_from values (phase 1)
# ---------------------------------------------------------------------------
#
# Each MAG has its own sample record in ENA. Its `related_sample_accession` field lists
# the samples the MAG was made from, for example:
#
#     SAMEA2622837:derived_from;SAMEA2622733:derived_from
#
# Submitters fill this in by hand, so the accessions come in four shapes:
#
#   SAMEA2622837            a BioSample, the usual case
#   ERS492926               the same sample under its older ENA accession. A search for one
#                           form does not find records that used the other, so we
#                           search with both.
#   ERR982800               a sequencing run instead of a sample. We look up which
#                           sample the run belongs to.
#   ERS9142775-ERS9142786   a range of samples. Kept exactly as written, not expanded,
#                           because we can't be sure the submitter used every sample
#                           in between.
#
# A MAG assembled from several samples lists all of them, in no fixed order.

PARENT_SAMPLE = "sample"
PARENT_RUN = "run"
PARENT_RANGE = "range"
PARENT_UNPARSABLE = "unparsable"

RELATION_DERIVED_FROM = "derived_from"

# What each kind of value in a derived_from field looks like:
#   RE_BIOSAMPLE   SAMEA2622837, SAMN39623922, SAMD00202456 (same pattern ENA uses)
#   RE_SECONDARY   ERS492926, also SRS... and DRS...
#   RE_RUN         ERR982800, also SRR... and DRR...
#   RE_RANGE       ERS9142775-ERS9142786
# The patterns are strict because these values are sent to ENA up to 500 per request,
# and a single malformed one (e.g. SAMN3962392, which is SAMN39623922 with a digit
# missing) makes ENA reject the whole request with HTTP 400.
RE_BIOSAMPLE = re.compile(r"^(SAMEA?[0-9]{6,}|SAM[ND][0-9]{8})$")
RE_SECONDARY = re.compile(r"^[EDSR]RS[0-9]{6,}$")
RE_RUN = re.compile(r"^[EDSR]RR[0-9]{6,}$")
RE_RANGE = re.compile(r"^[EDSR]RS[0-9]{6,}-[EDSR]RS[0-9]{6,}$")

# Spaces, tabs and quotes that submitters sometimes leave around a value,
# e.g. ' "SAMEA2622837" '.
_STRIP_CHARS = ' \t"\''


@dataclass(frozen=True)
class ParentValue:
    """One derived_from value as written (raw), and its shape (kind): PARENT_SAMPLE,
    PARENT_RUN, PARENT_RANGE or PARENT_UNPARSABLE."""

    raw: str
    kind: str


def classify_parent(token):
    """Return the shape of one derived_from value:

        SAMEA2622837, ERS492926   -> PARENT_SAMPLE
        ERR982800                 -> PARENT_RUN         its sample is looked up later
        ERS9142775-ERS9142786     -> PARENT_RANGE       kept as written, not expanded
        anything else             -> PARENT_UNPARSABLE
    """
    if RE_BIOSAMPLE.match(token) or RE_SECONDARY.match(token):
        return PARENT_SAMPLE
    if RE_RUN.match(token):
        return PARENT_RUN
    if RE_RANGE.match(token):
        return PARENT_RANGE
    return PARENT_UNPARSABLE


def parse_related_field(value):
    """Return the derived_from values of a related_sample_accession field.

        "SAMEA2622837:derived_from;SAMEA2622733:derived_from"
        -> [ParentValue("SAMEA2622837", PARENT_SAMPLE),
            ParentValue("SAMEA2622733", PARENT_SAMPLE)]

    Entries with any other relation are skipped. Values keep the order they were written
    in (fingerprint() sorts them later).

    A value of no known shape is kept as PARENT_UNPARSABLE, so the MAG can't look fully
    traced: it becomes partly_traced, or gets a status that says why.
    """
    if not value:
        return []
    out = []
    for chunk in value.split(";"):
        entry = chunk.strip(_STRIP_CHARS)
        if not entry or ":" not in entry:
            continue
        token, _, found_relation = entry.rpartition(":")
        token = token.strip(_STRIP_CHARS)
        if found_relation.strip(_STRIP_CHARS) != RELATION_DERIVED_FROM or not token:
            continue
        out.append(ParentValue(raw=token, kind=classify_parent(token)))
    return out


# ---------------------------------------------------------------------------
# Talking to ENA (phase 1)
# ---------------------------------------------------------------------------
#
# Every search is a POST to ENA_SEARCH_URL. A GET with hundreds of accessions makes the
# URL too long, and ENA rejects it.
#
# A request has failed when ENA returns an HTTP error, the connection fails, or the
# answer contains "ERROR occurred" (ENA sometimes reports a temporary failure this way,
# with HTTP 200). A failed request is retried. If it still fails, the run stops with
# an error, so a failure is never read as "ENA found nothing" and reported as
# "no duplicates".

ENA_SEARCH_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_ERROR_MARKER = "ERROR occurred"

# HTTP codes that mean the request itself is wrong, so it is not retried: 400 (e.g. an
# ERZ accession sent to the assembly endpoint, or SAMEA... and ERS... in one request),
# 404 (unknown endpoint), 414 (request too long). Other codes, such as 500 or 503, are
# retried.
PERMANENT_HTTP_CODES = frozenset({400, 404, 414})


class EnaError(RuntimeError):
    """A request ENA rejected, or one that still failed after all retries. The run stops
    with this error, so a failed request is never mistaken for "ENA found nothing"."""


class EnaHttpError(Exception):
    """An HTTP error status from ENA. EnaClient uses the code to decide whether to retry."""

    def __init__(self, status_code, body=""):
        super().__init__("HTTP %s" % status_code)
        self.status_code = status_code
        self.body = body


def _http_post(url, data, timeout):
    """POST form data and return the response text. An HTTP error status raises
    EnaHttpError."""
    request = urllib.request.Request(url, data=urllib.parse.urlencode(data).encode("utf-8"))
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as exc:
        try:
            body = exc.read().decode("utf-8", errors="replace")
        except OSError:
            body = ""
        raise EnaHttpError(exc.code, body)


def split_accessions_by_form(accessions):
    """Split accessions into the groups ENA accepts together in one request:

        SAMEA2622837, SAMN39623922    BioSample accessions
        ERS492926, SRS20266877        older sample accessions
        everything else               e.g. GCA_937896655, ERR982800

    A request with both SAMEA2622837 and ERS492926 in it gets HTTP 400.
    """
    biosample, secondary, other = [], [], []
    for accession in accessions:
        if RE_BIOSAMPLE.match(accession):
            biosample.append(accession)
        elif RE_SECONDARY.match(accession):
            secondary.append(accession)
        else:
            other.append(accession)
    return [group for group in (biosample, secondary, other) if group]


class EnaClient:
    """Sends searches to the ENA Portal API.

    A failed request is tried up to 4 times, waiting 5, 10 and 15 seconds in between.
    With cache_dir, each answer is saved as <cache_dir>/<hash of the request>.tsv, and
    the same request in a later run reads that file instead of asking ENA. request_count
    and cache_hits are shown in the phase 1 summary.
    """

    def __init__(self, cache_dir=None, timeout=180, max_retries=4, retry_wait=5.0,
                 batch_size=500):
        self.cache_dir = cache_dir
        self.timeout = timeout
        self.max_retries = max_retries
        self.retry_wait = retry_wait
        self.batch_size = batch_size
        self.request_count = 0
        self.cache_hits = 0
        if self.cache_dir:
            os.makedirs(self.cache_dir, exist_ok=True)

    def _cache_path(self, payload):
        """The cache file for this request, or None when there is no cache. Its name is the
        SHA-256 of the request, so the same request always maps to the same file."""
        if not self.cache_dir:
            return None
        key = json.dumps(payload, sort_keys=True)
        digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
        return os.path.join(self.cache_dir, digest + ".tsv")

    def _post(self, payload):
        """Send one request and return ENA's answer as text. Reads the cache first if there
        is one. Raises EnaError if the request is rejected or still fails after all
        retries."""
        path = self._cache_path(payload)
        if path and os.path.exists(path):
            self.cache_hits += 1
            with open(path, "r", encoding="utf-8") as handle:
                return handle.read()

        last_error = None
        for attempt in range(1, self.max_retries + 1):
            try:
                self.request_count += 1
                text = _http_post(ENA_SEARCH_URL, payload, self.timeout)
                if ENA_ERROR_MARKER in text:
                    last_error = "ENA's answer said '%s'" % ENA_ERROR_MARKER
                    LOG.warning("ENA's answer said '%s', attempt %s/%s (result=%s)",
                                ENA_ERROR_MARKER, attempt, self.max_retries,
                                payload.get("result"))
                else:
                    if path:
                        # Write to <file>.part first, then rename it. If the job is
                        # killed while writing, only the .part file is left, and the
                        # next run asks ENA again instead of reading a half-written answer.
                        partial = path + ".part"
                        with open(partial, "w", encoding="utf-8") as handle:
                            handle.write(text)
                        os.replace(partial, path)
                    return text
            except EnaHttpError as exc:
                last_error = "HTTP %s %s" % (exc.status_code, exc.body[:200])
                if exc.status_code in PERMANENT_HTTP_CODES:
                    raise EnaError(
                        "ENA rejected the request (result=%s): %s"
                        % (payload.get("result"), last_error)
                    )
                LOG.warning("ENA HTTP %s, attempt %s/%s", exc.status_code, attempt,
                            self.max_retries)
            except (OSError, http.client.HTTPException) as exc:
                # Network failures and timeouts. Anything else is a bug in this script and
                # is raised straight away instead of being retried.
                last_error = str(exc)
                LOG.warning("ENA request failed, attempt %s/%s: %s", attempt,
                            self.max_retries, exc)
            if attempt < self.max_retries and self.retry_wait:
                time.sleep(self.retry_wait * attempt)

        raise EnaError("ENA request failed after %s attempts (result=%s): %s"
                       % (self.max_retries, payload.get("result"), last_error))

    @staticmethod
    def _parse_tsv(text):
        """Turn ENA's TSV answer into one dict per row, keyed by the header line. A row with
        fewer values than the header is padded with "" so it is kept."""
        lines = [line for line in text.split("\n") if line.strip()]
        if len(lines) < 2:
            return []
        header = lines[0].split("\t")
        rows = []
        for line in lines[1:]:
            values = line.split("\t")
            values += [""] * (len(header) - len(values))
            rows.append(dict(zip(header, values)))
        return rows

    def search(self, result, fields, query=None, accessions=None):
        """Search one ENA result type (e.g. "assembly", "sample") and return every row as a
        dict with the requested fields.

        Give one of:
          query       a search expression, e.g. 'sample_accession="SAMEA2622837"'
          accessions  a list of accessions to fetch. Duplicates are removed, BioSamples
                      (SAM...) and older sample accessions (ERS...) go in separate
                      requests, and each request holds up to 500 accessions.
        """
        base = {
            "result": result,
            "fields": ",".join(fields),
            "format": "tsv",
            "limit": "0",  # 0 means no limit
        }

        if accessions is not None:
            rows = []
            for group in split_accessions_by_form(_dedupe(a for a in accessions if a)):
                for start in range(0, len(group), self.batch_size):
                    payload = dict(base)
                    payload["includeAccessions"] = ",".join(
                        group[start:start + self.batch_size])
                    rows.extend(self._parse_tsv(self._post(payload)))
            return rows

        payload = dict(base)
        payload["query"] = query or ""
        return self._parse_tsv(self._post(payload))


# ---------------------------------------------------------------------------
# Looking up each input genome in ENA (phase 1)
# ---------------------------------------------------------------------------

# Values of the source_sample_status column in mags.tsv:
#   traced                every derived_from link led to a source sample
#   partly_traced         some links led to source samples, but at least one did not, so
#                         the genome gets no source group
#   mag_sample_not_found  ENA lists no sample for the genome, or has no record of it
#                         (e.g. reference genomes such as GCA_000001405)
#   no_derived_from       the MAG sample's derived_from is empty or holds no accession
#                         (e.g. "missing")
#   source_not_found      derived_from names an accession, but no sample with runs was
#                         found behind it: ENA has no record of it (e.g. ERS950426), or
#                         it is a sample with no runs and nothing further back (e.g.
#                         SRS1615213)
#   trace_too_long        no sample with runs within MAX_TRACE_STEPS (10) steps
#   not_a_mag             ENA's assembly_type is in NON_MAG_ASSEMBLY_TYPES
#   repeated_input        the same genome as an earlier input line
#   accession_not_found   a GCA, ERZ or WGS accession ENA has no record of
#   invalid_accession     not a GCA, ERZ or WGS accession, e.g. GCF_000001405
STATUS_TRACED = "traced"
STATUS_PARTLY_TRACED = "partly_traced"
STATUS_MAG_SAMPLE_NOT_FOUND = "mag_sample_not_found"
STATUS_NO_DERIVED_FROM = "no_derived_from"
STATUS_SOURCE_NOT_FOUND = "source_not_found"
STATUS_TRACE_TOO_LONG = "trace_too_long"
STATUS_NOT_A_MAG = "not_a_mag"
STATUS_REPEATED_INPUT = "repeated_input"
STATUS_ACCESSION_NOT_FOUND = "accession_not_found"
STATUS_INVALID = "invalid_accession"

# ENA assembly_type values that are not MAGs. These genomes get not_a_mag and are not
# traced. Any other value, or none, counts as a MAG: a genome is left out only when ENA
# says it isn't one.
NON_MAG_ASSEMBLY_TYPES = frozenset({
    "clone or isolate",                                  # grown in culture
    "environmental single-cell amplified genome (sag)",  # from a single cell
    "primary metagenome",                                # a whole metagenome assembly
    "metatranscriptome",                                 # assembled from RNA
    "covid-19 outbreak",                                 # a SARS-CoV-2 genome
    "clinical isolate assembly",                         # an isolate from a patient
})

# The fields requested from each ENA endpoint. Each endpoint has different fields, and
# asking for one it doesn't have gets HTTP 400, so each has its own list.
ASSEMBLY_FIELDS = [
    "accession", "version", "study_accession", "center_name", "sample_accession",
    "run_accession", "binning_software", "assembly_software", "completeness_score",
    "contamination_score", "assembly_type", "last_updated", "tax_id", "scientific_name",
]
ANALYSIS_FIELDS = [
    "analysis_accession", "study_accession", "center_name",
    "sample_accession", "run_accession", "binning_software", "assembly_software",
    "completeness_score", "contamination_score", "assembly_type", "last_updated",
    "tax_id", "scientific_name",
    # Where ENA keeps an ERZ analysis's sequence: generated_ftp (e.g. for ERZ10779886)
    # or, when that is empty, submitted_ftp (e.g. for ERZ21930252). Used for the
    # download link in organism_check_candidates.tsv, because ENA's FASTA download
    # doesn't work for ERZ analyses.
    "generated_ftp", "submitted_ftp",
]
WGS_FIELDS = [
    "accession", "assembly_accession", "study_accession",
    "sample_accession", "run_accession", "assembly_type", "last_updated",
    "tax_id", "scientific_name",
]


@dataclass
class GenomeRecord:
    """One genome, as one row of mags.tsv. The field names are the column names.

    Read from ENA: the study, center, software, quality scores, taxonomy and fasta_ftp.
    Worked out by this script: sequence_check, sequence_md5, source_sample_status,
    source_samples, source_group_id, duplicate_status and gemsparcl_cluster (set in
    phase 2).

    version and assembly_accession are not written to mags.tsv: the version is added to
    genome_accession (GCA_937896655.1), and assembly_accession (the GCA of a WGS set,
    e.g. GCA_937896655 for CALBOS01) is only used to look the set up.
    """

    genome_accession: str
    input_accession: str
    source_sample_status: str
    version: str = ""
    mag_sample: str = ""
    study_accession: str = ""
    center_name: str = ""
    last_updated: str = ""
    assembly_software: str = ""
    binning_software: str = ""
    completeness_score_submitter: str = ""
    contamination_score_submitter: str = ""
    assembly_type: str = ""
    assembly_accession: str = ""
    tax_id: str = ""
    scientific_name: str = ""
    fasta_ftp: str = ""
    declared_runs: List[str] = field(default_factory=list)
    source_samples: List[str] = field(default_factory=list)
    source_group_id: str = ""
    duplicate_status: str = ""
    gemsparcl_cluster: str = ""
    sequence_check: str = ""
    sequence_md5: str = ""


# Where to look up each kind of input accession: (ENA endpoint, fields to request, column
# that holds the accession). For example, ERZ21930252 is searched in "analysis" and comes
# back in the "analysis_accession" column. An accession sent to the wrong endpoint gets
# HTTP 400.
ENDPOINTS = {
    KIND_GCA: ("assembly", ASSEMBLY_FIELDS, "accession"),
    KIND_ERZ: ("analysis", ANALYSIS_FIELDS, "analysis_accession"),
    KIND_WGS: ("wgs_set", WGS_FIELDS, "accession"),
}


# ENA fields copied into GenomeRecord as they are, under the same name.
_COPIED_FIELDS = [
    "version", "study_accession", "center_name", "last_updated", "assembly_software",
    "binning_software", "assembly_type", "assembly_accession", "tax_id", "scientific_name",
]
# The quality scores get a new name, to show in mags.tsv that they are the submitter's own.
_RENAMED_FIELDS = {
    "completeness_score": "completeness_score_submitter",
    "contamination_score": "contamination_score_submitter",
}


def _row_to_record(row, acc_key):
    """Turn one row from the assembly, analysis or wgs_set endpoint into a GenomeRecord.

    source_sample_status starts empty, and is set later to not_a_mag or to the result of
    tracing. run_accession can list several runs separated by ";" or ","; they become
    declared_runs.
    """
    declared = [r.strip() for r in re.split(r"[;,]", row.get("run_accession", ""))
                if r.strip()]
    return GenomeRecord(
        genome_accession=row.get(acc_key, ""),
        input_accession="",
        source_sample_status="",
        mag_sample=row.get("sample_accession", ""),
        declared_runs=declared,
        fasta_ftp=row.get("generated_ftp") or row.get("submitted_ftp", ""),
        **{name: row.get(name, "") for name in _COPIED_FIELDS},
        **{attr: row.get(name, "") for name, attr in _RENAMED_FIELDS.items()}
    )


def _version_key(value):
    """Turn a version into a number, so version 10 counts as newer than 9 (as text, "10"
    sorts before "9"). A missing version gives -1."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return -1


def _fetch_genomes(client, kind, accessions=None, query=None):
    """Fetch genomes of one kind (KIND_GCA, KIND_ERZ or KIND_WGS) by accession list or
    search query.

    ENA returns one row per version; only the latest is kept. Returns
    {unversioned accession: GenomeRecord}.
    """
    result, fields, acc_key = ENDPOINTS[kind]
    latest = {}
    for row in client.search(result, fields, query=query, accessions=accessions):
        record = _row_to_record(row, acc_key)
        current = latest.get(record.genome_accession)
        if current is None or _version_key(record.version) > _version_key(current.version):
            latest[record.genome_accession] = record
    return latest


def lookup_genomes(client, accessions):
    """Look up every input accession in ENA.

    Returns {input line: GenomeRecord}, one per line. A line that can't be looked up
    still gets a record, with source_sample_status saying why (invalid_accession,
    accession_not_found, repeated_input), so every input line appears in mags.tsv.
    """
    records = {}
    wanted = []
    for accession in accessions:
        if accession.kind == KIND_INVALID:
            records[accession.raw] = GenomeRecord(
                genome_accession=accession.raw,
                input_accession=accession.raw,
                source_sample_status=STATUS_INVALID,
            )
        else:
            wanted.append(accession)

    if not wanted:
        return records

    # Look up each distinct accession once: {accession sent to ENA: record}.
    found = {}
    for kind in ENDPOINTS:
        queries = sorted({a.query for a in wanted if a.kind == kind})
        if queries:
            found.update(_fetch_genomes(client, kind, accessions=queries))

    # The wgs_set endpoint has no center, software or quality fields, so a WGS set is
    # replaced by the record of its GCA (its assembly_accession).
    wgs_gcas = {a.query: found[a.query].assembly_accession.strip()
                for a in wanted
                if a.kind == KIND_WGS and a.query in found
                and found[a.query].assembly_accession.strip()}
    if wgs_gcas:
        upgraded = _fetch_genomes(client, KIND_GCA, accessions=sorted(set(wgs_gcas.values())))
        for query, gca in wgs_gcas.items():
            if gca in upgraded:
                found[query] = upgraded[gca]

    # One record per input line. When two lines give the same genome (GCA_937896655 and
    # GCA_937896655.1, or CALBOS01 and its GCA_937896655), the later line gets
    # repeated_input.
    taken = set()
    for accession in wanted:
        record = found.get(accession.query)
        if record is None:
            record = GenomeRecord(genome_accession=accession.query, input_accession="",
                                  source_sample_status=STATUS_ACCESSION_NOT_FOUND)
        if record.genome_accession in taken:
            record = GenomeRecord(genome_accession=record.genome_accession,
                                  version=record.version, input_accession="",
                                  source_sample_status=STATUS_REPEATED_INPUT)
        taken.add(record.genome_accession)
        records[accession.raw] = replace(record, input_accession=accession.raw)

    _mark_non_mags(records.values())
    return records


def _mark_non_mags(records):
    """Set source_sample_status to not_a_mag for genomes whose assembly_type is in
    NON_MAG_ASSEMBLY_TYPES."""
    for record in records:
        if (not record.source_sample_status
                and (record.assembly_type or "").lower() in NON_MAG_ASSEMBLY_TYPES):
            record.source_sample_status = STATUS_NOT_A_MAG


# ---------------------------------------------------------------------------
# Identical sequences (phase 1, step 1, with --fasta-dir)
# ---------------------------------------------------------------------------
#
# Genomes with the same sequence are compared only once: one copy is kept for the source
# groups and gemsparcl, and the others get duplicate_status duplicate_sequence. Every
# copy is still traced, so its provenance is recorded.
#
# A genome's checksum is the MD5 of the sorted MD5s of its sequences (uppercased, without
# whitespace), so headers, contig order, line length and case don't change it.

# Values of the sequence_check column (empty if not checked: no --fasta-dir, an
# invalid, not found, not_a_mag or repeated input, or a genome found by reverse lookup):
#   unique          no other checked genome has the same sequence
#   duplicate_kept  the copy kept for the comparisons (see identical_sequences.tsv)
#   duplicate       another copy; not compared again
#   no_fasta_file   no usable FASTA file: missing, another version, empty or unreadable
# The sequence_md5 column holds each checked genome's checksum, for matching in later
# runs.
SEQ_UNIQUE = "unique"
SEQ_DUPLICATE_KEPT = "duplicate_kept"
SEQ_DUPLICATE = "duplicate"
SEQ_NO_FASTA_FILE = "no_fasta_file"

FASTA_EXTENSIONS = (".fa", ".fasta", ".fna", ".fsa")
IDENTICAL_FILE = "identical_sequences.tsv"
IDENTICAL_COLUMNS = ["sequence_md5", "kept_genome", "n_genomes", "genomes", "studies"]

# A copy from EMG is preferred as the kept one.
EBI_CENTER_NAME = "EMG"


def _accession_of_file(name):
    """Read the accession from a FASTA file name:

        GCA_937896655.1.fna.gz  -> ("GCA_937896655", "1")
        ERZ21930252.fa          -> ("ERZ21930252", "")
        notes.txt               -> None
    """
    if name.endswith(".gz"):
        name = name[:-3]
    for extension in FASTA_EXTENSIONS:
        if name.endswith(extension):
            return split_version(name[:-len(extension)])
    return None


def find_fasta_files(fasta_dir, records):
    """Find each record's FASTA file in fasta_dir by its name: GCA_937896655.1.fna.gz or
    GCA_937896655.fna.gz for GCA_937896655 at version 1. A file named with another
    version, such as GCA_937896655.2.fna.gz, is not used. Returns [(record, path)] and
    the names of the files not used."""
    files_by_accession = {}  # {accession: [(version, file name)]}
    for name in sorted(os.listdir(fasta_dir)):
        parsed = _accession_of_file(name)
        if parsed:
            files_by_accession.setdefault(parsed[0], []).append((parsed[1], name))
    found, other_version = [], []
    for record in records:
        files = files_by_accession.get(record.genome_accession, [])
        usable = [name for version, name in files
                  if not version or not record.version or version == record.version]
        if usable:
            found.append((record, os.path.join(fasta_dir, usable[0])))
        else:
            other_version.extend(name for _, name in files)
    return found, other_version


def _read_fasta_sequences(path):
    """Read the sequences of a FASTA file (plain or gzipped) one at a time, uppercased
    and without whitespace."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        parts = None
        for line in handle:
            if line.startswith(">"):
                if parts is not None:
                    yield "".join(parts).upper()
                parts = []
            elif parts is not None:
                parts.append("".join(line.split()))
        if parts is not None:
            yield "".join(parts).upper()


def sequence_checksum(path):
    """Return the genome's checksum: the MD5 of the sorted MD5s of its sequences, or ""
    if the file holds no sequence."""
    hashes = sorted(hashlib.md5(sequence.encode("ascii", errors="replace")).hexdigest()
                    for sequence in _read_fasta_sequences(path) if sequence)
    if not hashes:
        return ""
    return hashlib.md5("\n".join(hashes).encode("ascii")).hexdigest()


def find_identical_sequences(records, fasta_dir):
    """Record each input MAG's checksum in sequence_md5, and set sequence_check to unique
    or no_fasta_file. Returns the sets of identical genomes as {checksum: genomes};
    choose_kept_copies() picks the kept copy after tracing."""
    files, other_version = find_fasta_files(fasta_dir, records)
    if other_version:
        LOG.warning("%s FASTA file(s) in %s hold another version of the genome than ENA's "
                    "current one and were not used, e.g. %s",
                    _number(len(other_version)), fasta_dir, other_version[0])
    for record in records:
        record.sequence_check = SEQ_NO_FASTA_FILE
    by_checksum = {}
    for record, path in files:
        try:
            checksum = sequence_checksum(path)
        except (OSError, EOFError) as exc:
            LOG.warning("can't read %s: %s", path, exc)
            continue
        if not checksum:
            LOG.warning("%s holds no sequence", path)
            continue
        record.sequence_md5 = checksum
        by_checksum.setdefault(checksum, []).append(record)
    for same in by_checksum.values():
        if len(same) == 1:
            same[0].sequence_check = SEQ_UNIQUE
    return {checksum: same for checksum, same in by_checksum.items() if len(same) >= 2}


def _kept_first(record):
    """The order in which copies are chosen as the kept one: a traced copy first (only
    traced genomes get a source group), then a copy from EMG, then a GCA before an ERZ
    (the GCA is the INSDC accession, with a version), then the lowest accession."""
    return (record.source_sample_status != STATUS_TRACED,
            record.center_name != EBI_CENTER_NAME,
            not record.genome_accession.startswith("GCA_"),
            record.genome_accession)


def choose_kept_copies(sets):
    """After tracing, keep one copy of each set of identical genomes (see _kept_first) and
    mark the others duplicate. Each set is sorted with the kept copy first."""
    for same in sets.values():
        same.sort(key=_kept_first)
        same[0].sequence_check = SEQ_DUPLICATE_KEPT
        for record in same[1:]:
            record.sequence_check = SEQ_DUPLICATE


def write_identical_sequences_tsv(path, sets):
    """Write identical_sequences.tsv: one row per set of identical genomes, with the
    kept genome listed first."""
    rows = []
    for checksum, same in sorted(sets.items(), key=lambda item: versioned(item[1][0])):
        rows.append([checksum, versioned(same[0]), str(len(same)),
                     _join(versioned(r) for r in same),
                     _join(_dedupe(r.study_accession for r in same if r.study_accession))])
    _write_tsv(path, IDENTICAL_COLUMNS, rows)


# ---------------------------------------------------------------------------
# Tracing each MAG sample back to its source samples (phase 1)
# ---------------------------------------------------------------------------
#
# A MAG has its own sample record in ENA, but the sequencing data is not stored under it.
# Its derived_from field names the sample(s) the MAG was built from. The script follows
# these links until it reaches a sample that has runs (the sequencing data): that sample
# is a source sample. A link can also name a run directly; its sample is then the source
# sample.
#
#   SAMEA112717052 (MAG sample) -> SAMEA3504092 (has run ERR982800)
#   SAMEA112374585 (MAG sample) -> ERR982800 (a run of SAMEA3504092)
#   Both MAGs have the source sample SAMEA3504092.
#
# A linked sample without runs is followed further back through its own derived_from.
# A sample with runs ends the trace, even if it has a derived_from of its own, because
# that link points to where the material came from (e.g. the host), not to more data.
# When derived_from names several samples, each is followed, and the MAG gets all the
# source samples they lead to.
#
# Each step looks up all MAG samples together (one search for runs, one for samples),
# not one MAG at a time.

SAMPLE_FIELDS = ["sample_accession", "secondary_sample_accession",
                 "related_sample_accession"]
READ_RUN_FIELDS = ["run_accession", "sample_accession", "secondary_sample_accession"]

# The most derived_from steps followed back from a MAG sample, so a very long chain can't
# keep sending requests. Loops (A's derived_from names B, and B's names A) end without
# this limit: a sample already looked up for a MAG is not looked up again.
MAX_TRACE_STEPS = 10


@dataclass
class TraceResult:
    """What was found by tracing one MAG sample, e.g. SAMEA112717052:

    source_samples  ["SAMEA3504092"]
    search_keys     {"SAMEA3504092", "ERS811241"}: every accession of the source samples
    runs            ["ERR982800"]: the runs of the source samples
    status          "traced": the source_sample_status

    search_keys and runs are used only for the reverse lookup: finding other genomes
    from the same source samples.
    """

    source_samples: List[str] = field(default_factory=list)
    search_keys: Set[str] = field(default_factory=set)
    runs: List[str] = field(default_factory=list)
    status: str = ""


def _both_forms(row):
    """The two accessions ENA gives for one sample: BioSample, and INSDC secondary."""
    return [row.get("sample_accession", "").strip(),
            row.get("secondary_sample_accession", "").strip()]


def _canonical_sample(row, fallback=""):
    """Return the one accession we use for this sample, so the same sample is always
    written the same way: its BioSample (SAMEA2622837) if ENA gives one, otherwise its
    older accession (ERS492926), otherwise fallback (the value as written in derived_from).
    """
    biosample, secondary = _both_forms(row)
    return biosample or secondary or fallback


def _index_samples_by_either_form(rows):
    """Index sample rows under both of their accessions, so a lookup by either one hits."""
    index = {}
    for row in rows:
        for form in _both_forms(row):
            if form:
                index.setdefault(form, row)
    return index


def _index_runs(rows):
    """Make two lookups from read_run rows: by run (ERR982800 -> its row) and by sample
    (SAMEA3504092 or ERS811241 -> the rows of all its runs)."""
    by_run = {}
    by_sample = {}
    for row in rows:
        run = row.get("run_accession", "").strip()
        if run:
            by_run.setdefault(run, row)
        for form in _both_forms(row):
            if form:
                by_sample.setdefault(form, []).append(row)
    return by_run, by_sample


class _Trace:
    """The tracing of one MAG sample: what it has found so far (result), the links to
    look up in the next step (to_look_up), and the samples already looked up (visited)."""

    def __init__(self, mag_sample):
        self.result = TraceResult()
        self.visited = {mag_sample}
        self.to_look_up = []
        # A derived_from value that is not an accession (e.g. "missing").
        self.saw_unparsable = False
        # A link that leads to no data: a run or sample ENA has no record of, or a sample
        # with no runs and nothing further back.
        self.saw_missing_parent = False
        # ENA has no record of the MAG sample itself.
        self.mag_sample_missing = False

    def add_source(self, canonical, search_keys=(), runs=()):
        if canonical and canonical not in self.result.source_samples:
            self.result.source_samples.append(canonical)
        self.result.search_keys.update(form for form in search_keys if form)
        for run in runs:
            if run and run not in self.result.runs:
                self.result.runs.append(run)

    def enqueue(self, parent):
        """Decide what to do with one derived_from value of this MAG's sample:

            ERS9142775-ERS9142786     a range: recorded as a source sample, as written
            missing                   not an accession: noted, so the status can say why
                                      tracing failed
            SAMEA3504092, ERR982800   a sample or run: looked up in the next step,
                                      unless already looked up for this MAG
        """
        if parent.kind == PARENT_RANGE:
            # Ranges are not expanded, so only MAGs that list exactly the same range share
            # a source group.
            self.add_source(parent.raw, search_keys=[parent.raw])
        elif parent.kind == PARENT_UNPARSABLE:
            self.saw_unparsable = True
        elif parent.raw not in self.visited:
            self.visited.add(parent.raw)
            self.to_look_up.append(parent)

    def finish(self):
        """Set the status and return the result.

        If at least one source sample was found, the status is traced, or partly_traced
        if some other link could not be followed.

        If none was found, the status gives the reason. When more than one reason
        applies, the first of these is used: mag_sample_not_found, trace_too_long,
        source_not_found, no_derived_from.
        """
        unfinished = bool(self.to_look_up)
        if self.result.source_samples:
            partial = self.saw_unparsable or self.saw_missing_parent or unfinished
            self.result.status = STATUS_PARTLY_TRACED if partial else STATUS_TRACED
        elif self.mag_sample_missing:
            self.result.status = STATUS_MAG_SAMPLE_NOT_FOUND
        elif unfinished:
            self.result.status = STATUS_TRACE_TOO_LONG
        elif self.saw_missing_parent:
            self.result.status = STATUS_SOURCE_NOT_FOUND
        else:
            self.result.status = STATUS_NO_DERIVED_FROM
        return self.result


def trace_source_samples(client, mag_samples):
    """Trace each MAG sample back to its source samples.

    Returns {MAG sample: TraceResult}, with an entry for every MAG sample, including
    those that could not be traced (their status says why).
    """
    wanted = _dedupe(s.strip() for s in mag_samples if s and s.strip())
    if not wanted:
        return {}

    traces = {sample: _Trace(sample) for sample in wanted}

    # Fetch the ENA records of all MAG samples, up to 500 per request, and read their
    # derived_from links.
    index = _index_samples_by_either_form(
        client.search("sample", SAMPLE_FIELDS, accessions=wanted))
    for mag_sample, trace in traces.items():
        row = index.get(mag_sample)
        if row is None:
            trace.mag_sample_missing = True
            continue
        trace.visited.update(form for form in _both_forms(row) if form)
        for parent in parse_related_field(row.get("related_sample_accession")):
            trace.enqueue(parent)

    for depth in range(1, MAX_TRACE_STEPS + 1):
        step_links = {parent.raw for trace in traces.values()
                      for parent in trace.to_look_up}
        if not step_links:
            break
        LOG.debug("tracing source samples, step %s: looking up %s accession(s)",
                  depth, len(step_links))

        by_run, by_sample = _index_runs(
            client.search("read_run", READ_RUN_FIELDS, accessions=sorted(step_links)))

        without_runs = {}
        for mag_sample, trace in traces.items():
            current, trace.to_look_up = trace.to_look_up, []
            for parent in current:
                if parent.kind == PARENT_RUN:
                    row = by_run.get(parent.raw)
                    if row is None or not _canonical_sample(row):
                        # The link names a run ENA has no record of. A run has no
                        # derived_from of its own, so this link can't be followed further.
                        trace.saw_missing_parent = True
                        continue
                    # Look up the run's sample in the next step. That also gives the
                    # sample's other runs, which other MAGs may name instead: e.g.
                    # SAMEA112376089 names ERR1614694, one of SAMEA4427029's 8 runs.
                    trace.enqueue(ParentValue(_canonical_sample(row), PARENT_SAMPLE))
                    continue
                rows = by_sample.get(parent.raw)
                if rows:
                    trace.add_source(
                        _canonical_sample(rows[0], fallback=parent.raw),
                        search_keys=_both_forms(rows[0]) + [parent.raw],
                        runs=[r.get("run_accession", "").strip() for r in rows],
                    )
                else:
                    without_runs.setdefault(parent.raw, []).append(mag_sample)

        if not without_runs:
            continue
        # Samples without runs: read their own derived_from links, to follow in the next step.
        parents = _index_samples_by_either_form(
            client.search("sample", SAMPLE_FIELDS, accessions=sorted(without_runs)))
        for token in sorted(without_runs):
            row = parents.get(token)
            for mag_sample in without_runs[token]:
                trace = traces[mag_sample]
                if row is None:
                    trace.saw_missing_parent = True
                    continue
                trace.visited.update(form for form in _both_forms(row) if form)
                before = (len(trace.to_look_up), len(trace.result.source_samples))
                for parent in parse_related_field(row.get("related_sample_accession")):
                    trace.enqueue(parent)
                if (len(trace.to_look_up), len(trace.result.source_samples)) == before:
                    # A sample with no runs and nothing new to follow (e.g. SRS1615213):
                    # this link leads to no data.
                    trace.saw_missing_parent = True

    unfinished = sum(1 for trace in traces.values() if trace.to_look_up)
    if unfinished:
        LOG.warning("%s MAG sample(s) still had links to follow after %s steps; they are "
                    "marked trace_too_long or partly_traced",
                    unfinished, MAX_TRACE_STEPS)

    return {mag_sample: trace.finish() for mag_sample, trace in traces.items()}


# ---------------------------------------------------------------------------
# Reverse lookup: finding other genomes from the same source samples (phase 1)
# ---------------------------------------------------------------------------
#
# Tracing goes back from a MAG to its source samples. The reverse lookup goes the other
# way, from the source samples to every genome built from them, because the genomes
# other studies built from them are usually not in the input list. ENA is searched for
# every sample whose derived_from names a source sample. ENA matches only the exact
# accession written, and submitters name the same sample in different ways, so each of
# them is searched. For SAMEA2622837:
#
#   SAMEA2622837   its BioSample accession
#   ERS492926      its older accession
#   ERR599156      one of its runs
#
# A source sample that is a range (e.g. ERS9142775-ERS9142786) is searched as written.
#
# A sample found this way that has no MAG (e.g. a bin) becomes a search key in the next
# step, because MAGs can name the bin instead of the source sample. This repeats for up
# to MAX_TRACE_STEPS steps.

def build_derived_from_query(keys):
    """Build the ENA search for samples whose derived_from names any of keys:

        ["SAMEA2622837", "ERS492926"] ->
        related_sample_accession="SAMEA2622837:derived_from" OR
        related_sample_accession="ERS492926:derived_from"

    Unlike a list of accessions, a search can mix SAMEA and ERS accessions."""
    terms = ['related_sample_accession="%s:derived_from"' % key for key in keys]
    return " OR ".join(terms)


def source_search_keys(trace_results):
    """The accessions and runs of all source samples found, for the reverse lookup. This
    includes partly_traced genomes: some of their links failed, but the source samples
    they have are correct."""
    keys = set()
    for trace_result in trace_results:
        if trace_result.status not in (STATUS_TRACED, STATUS_PARTLY_TRACED):
            continue
        keys.update(form for form in trace_result.search_keys if form)
        keys.update(run for run in trace_result.runs if run)
    return sorted(keys)


def find_samples_derived_from(client, keys):
    """Search ENA for every sample whose derived_from names one of keys. Returns
    {sample: its accessions}, in the order found; each sample is named by its BioSample
    accession where ENA has one."""
    found = {}
    for start in range(0, len(keys), client.batch_size):
        query = build_derived_from_query(keys[start:start + client.batch_size])
        for row in client.search("sample", SAMPLE_FIELDS, query=query):
            accession = _canonical_sample(row)
            if accession:
                found.setdefault(accession, [form for form in _both_forms(row) if form])
    return found


def find_other_genomes(client, keys, known_samples):
    """The reverse lookup: search ENA for genomes built from the source samples, step by
    step as described above. keys are the source samples' accessions and runs. The input
    genomes' own samples (known_samples) are skipped."""
    seen = set(known_samples)
    genomes = []
    for depth in range(1, MAX_TRACE_STEPS + 1):
        if not keys:
            break
        found = {sample: accessions
                 for sample, accessions in find_samples_derived_from(client, keys).items()
                 if sample not in seen}
        LOG.debug("reverse lookup, step %s: %s accession(s) searched, %s new "
                  "sample(s) found", depth, len(keys), len(found))
        seen.update(found)
        step_genomes, without_mag = find_genomes_of_samples(client, list(found))
        genomes.extend(step_genomes)
        keys = sorted(accession for sample in without_mag for accession in found[sample])
    return genomes


def find_genomes_of_samples(client, samples):
    """Find the genomes of samples found by find_samples_derived_from(). Both the
    assembly (GCA) and analysis (ERZ) endpoints are searched, because another study's
    genome can be either. Returns the genomes (non-MAGs included, as not_a_mag) and the
    samples with no MAG."""
    by_sample = {}
    for start in range(0, len(samples), client.batch_size):
        chunk = samples[start:start + client.batch_size]
        query = " OR ".join('sample_accession="%s"' % s for s in chunk)
        for kind in (KIND_GCA, KIND_ERZ):
            for record in _fetch_genomes(client, kind, query=query).values():
                by_sample.setdefault(record.mag_sample, []).append(record)

    genomes, without_mag = [], []
    for sample in samples:
        records = by_sample.get(sample, [])
        _mark_non_mags(records)
        mags = [r for r in records if r.source_sample_status != STATUS_NOT_A_MAG]
        # A MAG submitted as an ERZ analysis also has a GCA under the same sample (e.g.
        # ERZ21930820 and GCA_963687975 under SAMEA114724776), so ERZs are kept only
        # when the sample has no GCA.
        gcas = [r for r in mags if r.genome_accession.startswith("GCA_")]
        genomes.extend(gcas or mags)
        genomes.extend(r for r in records if r.source_sample_status == STATUS_NOT_A_MAG)
        if not mags:
            without_mag.append(sample)
    return genomes, without_mag


# ---------------------------------------------------------------------------
# Source groups (phase 1)
# ---------------------------------------------------------------------------
#
# Genomes with exactly the same source samples form a source group. Only groups with
# genomes from 2 or more studies can hold duplicates; phase 2 checks those.

# Values of the duplicate_status column. Phase 1 sets to_compare, single_study,
# no_source_group or duplicate_sequence. Phase 2 looks again at the to_compare genomes
# and sets duplicate, not_duplicate, other_study_no_cluster or no_cluster.
#   duplicate               a genome from another study in its source group has the
#                           same gemsparcl cluster
#   not_duplicate           the genomes from other studies all have other clusters
#   other_study_no_cluster  no genome from another study has its cluster, but one has no
#                           cluster, so this could still be a duplicate
#   no_cluster              this genome is not in the gemsparcl clusters file
#   to_compare              its source group has genomes from 2 or more studies, so it
#                           is compared in phase 2
#   single_study            all genomes in its source group are from one study, so it
#                           can't be a duplicate
#   no_source_group         its source_sample_status is not traced, so it can't be
#                           compared
#   duplicate_sequence      has the same sequence as the duplicate_kept genome, which is
#                           compared instead
DUP_DUPLICATE = "duplicate"
DUP_NOT_DUPLICATE = "not_duplicate"
DUP_OTHER_STUDY_NO_CLUSTER = "other_study_no_cluster"
DUP_NO_CLUSTER = "no_cluster"
DUP_TO_COMPARE = "to_compare"
DUP_SINGLE_STUDY = "single_study"
DUP_NO_SOURCE_GROUP = "no_source_group"
DUP_DUPLICATE_SEQUENCE = "duplicate_sequence"


def _short_id(prefix, items):
    """A short id that stays the same in every run: prefix plus 12 characters of a hash
    of items."""
    return prefix + hashlib.sha1("|".join(items).encode("utf-8")).hexdigest()[:12]


def fingerprint(source_samples):
    """Return the source group id, made from the sorted source samples, so the same
    samples always give the same id:

        ["SAMEA2622837", "SAMEA2622733"] -> sg_3bd7455b797d

    Sorting matters because the same samples can be listed in either order, even within
    one study: in PRJEB82487, SAMEA117322715 lists SAMEA2622837;SAMEA2622733 and
    SAMEA117322723 lists SAMEA2622733;SAMEA2622837. Both get sg_3bd7455b797d."""
    return _short_id("sg_", sorted(set(sample for sample in source_samples if sample)))


def build_source_groups(records):
    """Give every fully traced genome its source group id and its phase 1
    duplicate_status. Identical copies (sequence_check duplicate) get no group, so they
    are not compared again. Returns {source_group_id: set of study accessions}."""
    for record in records:
        if record.sequence_check == SEQ_DUPLICATE:
            record.source_group_id = ""
            record.duplicate_status = DUP_DUPLICATE_SEQUENCE
        elif record.source_sample_status == STATUS_TRACED and record.source_samples:
            record.source_group_id = fingerprint(record.source_samples)
        else:
            record.source_group_id = ""
            record.duplicate_status = DUP_NO_SOURCE_GROUP

    studies_by_group_id = studies_by_group(records)
    for record in records:
        if record.source_group_id:
            multi = _multi_study(studies_by_group_id, record.source_group_id)
            record.duplicate_status = DUP_TO_COMPARE if multi else DUP_SINGLE_STUDY
    return studies_by_group_id


def studies_by_group(records):
    """Return {source_group_id: set of study accessions} for records that have a group id.

    Phase 2 uses this to rebuild the groups from mags.tsv.
    """
    studies = {}
    for record in records:
        if not record.source_group_id:
            continue
        group_studies = studies.setdefault(record.source_group_id, set())
        if record.study_accession:
            group_studies.add(record.study_accession)
    return studies


def _multi_study(studies_by_group_id, group_id):
    return len(studies_by_group_id.get(group_id, ())) >= 2


# ---------------------------------------------------------------------------
# Output (phases 1 and 2)
# ---------------------------------------------------------------------------

MAGS_COLUMNS = [
    "source_group_id", "genome_accession", "input_accession",
    "mag_sample", "source_samples",
    "declared_runs",
    "study_accession", "center_name", "last_updated",
    "assembly_software", "binning_software",
    "completeness_score_submitter", "contamination_score_submitter",
    "assembly_type",
    "tax_id", "scientific_name", "fasta_ftp",
    "gemsparcl_cluster",
    "duplicate_status", "source_sample_status", "sequence_check", "sequence_md5",
]


def _join(values):
    return ",".join(values)


def versioned(record):
    """GCA_937896655 at version 1 is written as GCA_937896655.1."""
    if record.version:
        return record.genome_accession + "." + record.version
    return record.genome_accession


def split_version(accession):
    """GCA_937896655.1 -> ("GCA_937896655", "1")."""
    base, dot, version = accession.rpartition(".")
    if dot and version.isdigit():
        return base, version
    return accession, ""


def _to_cell(value):
    """A list is written sorted and joined by commas:
    ["SAMEA2622837", "SAMEA2622733"] -> "SAMEA2622733,SAMEA2622837"."""
    if isinstance(value, list):
        return _join(sorted(set(value)))
    return value


def _from_cell(field_type, text):
    """The reverse of _to_cell, for list fields; other fields stay as text."""
    if field_type == List[str]:
        return [part for part in text.split(",") if part]
    return text


def _write_tsv(path, columns, rows):
    """Write a tab-separated file: a header line with the column names, then the rows."""
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def write_mags_tsv(path, records):
    """Write one row per genome, including inputs that could not be traced (their
    source_sample_status says why), so every input line appears in the output."""
    rows = []
    for record in sorted(records, key=lambda r: (r.source_group_id, r.genome_accession,
                                                 r.input_accession)):
        cells = {name: _to_cell(getattr(record, name)) for name in MAGS_COLUMNS}
        cells["genome_accession"] = versioned(record)
        rows.append([cells[name] for name in MAGS_COLUMNS])
    _write_tsv(path, MAGS_COLUMNS, rows)


def read_mags_tsv(path):
    """Rebuild GenomeRecords from mags.tsv. A column with no GenomeRecord attribute of
    the same name (e.g. one added by hand) is skipped."""
    types = {f.name: f.type for f in dataclass_fields(GenomeRecord)}
    records = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            record = GenomeRecord(**{name: _from_cell(types[name], text or "")
                                     for name, text in row.items() if name in types})
            record.genome_accession, record.version = split_version(record.genome_accession)
            records.append(record)
    return records


def _number(value):
    """1234 -> "1,234"."""
    return "{:,}".format(value)


def _log_row(label, count, note=""):
    """One line of the run summary: a label, a right-aligned count and an optional note."""
    LOG.info("  %-36s %7s%s", label, _number(count), "  (%s)" % note if note else "")


def _log_status_table(title, records, statuses):
    """Count the records with each duplicate_status in statuses, in two columns: genomes
    from the input list, and genomes found by reverse lookup."""
    inputs = Counter(r.duplicate_status for r in records if r.input_accession)
    found = Counter(r.duplicate_status for r in records if not r.input_accession)
    table = "  %-36s %14s %16s"
    LOG.info("")
    LOG.info(table, title, "in your list", "reverse lookup")
    for status in statuses:
        LOG.info(table, "  " + status, _number(inputs[status]), _number(found[status]))
    LOG.info(table, "  total", _number(sum(inputs[s] for s in statuses)),
             _number(sum(found[s] for s in statuses)))


# Inputs that can't be traced, so they are left out of the traced rate.
UNTRACEABLE_STATUSES = (STATUS_INVALID, STATUS_ACCESSION_NOT_FOUND, STATUS_NOT_A_MAG,
                        STATUS_REPEATED_INPUT)


def _log_sequence_check(inputs, identical):
    """Print step 1 of the phase 1 summary. identical holds the sets of identical genomes,
    or is None without --fasta-dir."""
    LOG.info("Step 1: compare the input genomes' sequences")
    if identical is None:
        LOG.info("  not run (no --fasta-dir)")
        return
    counts = Counter(r.sequence_check for r in inputs)
    LOG.info("  sequence_check of the inputs:")
    for value in (SEQ_UNIQUE, SEQ_DUPLICATE_KEPT, SEQ_DUPLICATE, SEQ_NO_FASTA_FILE):
        note = ""
        if value == SEQ_DUPLICATE_KEPT:
            note = "one per set of identical sequences, in %s" % IDENTICAL_FILE
        _log_row("  " + value, counts[value], note)


def log_phase1_summary(records, studies_by_group_id, client, identical=None):
    """Print the phase 1 summary. identical holds the sets of identical genomes, or is
    None without --fasta-dir."""
    inputs = [r for r in records if r.input_accession]
    statuses = Counter(r.source_sample_status for r in inputs)
    traceable = sum(n for s, n in statuses.items() if s not in UNTRACEABLE_STATUSES)
    left_out = ["%s %s" % (_number(n), s) for s, n in statuses.most_common()
                if s in UNTRACEABLE_STATUSES]
    rate = "%.1f%% of the inputs" % (
        100.0 * statuses[STATUS_TRACED] / traceable if traceable else 0.0)
    if left_out:
        rate += ", leaving out " + ", ".join(left_out)
    multi = sum(1 for studies in studies_by_group_id.values() if len(studies) >= 2)

    _log_row("input accessions", len(inputs))
    _log_sequence_check(inputs, identical)
    LOG.info("Step 2: trace genomes to their source samples, then find every genome "
             "built from those samples (reverse lookup)")
    LOG.info("  source_sample_status of the inputs:")
    for status, count in statuses.most_common():
        _log_row("  " + status, count, rate if status == STATUS_TRACED else "")
    _log_row("genomes found by reverse lookup", len(records) - len(inputs),
             "from the same source samples, not in your list")
    _log_row("source groups", len(studies_by_group_id),
             "%s used by 2 or more studies" % _number(multi))
    _log_row("ENA requests", client.request_count,
             "%s answered from the cache" % _number(client.cache_hits))
    _log_status_table("duplicate_status after phase 1", records,
                      (DUP_TO_COMPARE, DUP_SINGLE_STUDY, DUP_DUPLICATE_SEQUENCE,
                       DUP_NO_SOURCE_GROUP))


# ---------------------------------------------------------------------------
# Reading the gemsparcl clusters (phase 2)
# ---------------------------------------------------------------------------
#
# gemsparcl assigns genomes of the same organism to the same cluster. The script reads
# the *_clusters.csv of a `gemsparcl cluster` run; it never runs gemsparcl itself.

ENA_FASTA_URL = "https://www.ebi.ac.uk/ena/browser/api/fasta/%s?download=true"


def read_gemsparcl_clusters(path, genomes):
    """Read a gemsparcl *_clusters.csv into {genome accession: cluster}, keeping only the
    given genomes (versioned, e.g. GCA_937896655.1), because the file can hold millions
    of rows. An id with a version matches only that version.

    An id without a version (GCA_937896655) matches the genome's version in mags.tsv.
    This is only for benchmarking with older clusters files that have no versions, and
    the number of genomes matched this way is logged."""
    by_unversioned = {split_version(genome)[0]: genome for genome in genomes}
    assignments = {}
    matched_without_version = set()
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        columns = {name.lower(): name for name in (reader.fieldnames or ())}
        id_col = columns.get("genome_id")
        cluster_col = columns.get("cluster")
        if not id_col or not cluster_col:
            raise ValueError(
                "%s: expected genome_id and cluster columns, found %s"
                % (path, reader.fieldnames))
        for row in reader:
            genome = (row.get(id_col) or "").strip()
            cluster = (row.get(cluster_col) or "").strip()
            if not genome or not cluster:
                continue
            if genome in genomes:
                assignments[genome] = cluster
                matched_without_version.discard(genome)
            elif genome in by_unversioned and by_unversioned[genome] not in assignments:
                assignments[by_unversioned[genome]] = cluster
                matched_without_version.add(by_unversioned[genome])
    if matched_without_version:
        LOG.warning("%s genome(s) were matched without a version: %s has ids such as "
                    "GCA_937896655 instead of GCA_937896655.1",
                    _number(len(matched_without_version)), path)
    return assignments


def attach_clusters(records, assignments):
    """Set each genome's gemsparcl_cluster from the clusters file."""
    for record in records:
        record.gemsparcl_cluster = assignments.get(versioned(record), "")


def _fasta_url(record):
    """The link to download the genome's sequence: fasta_ftp if ENA has one (ERZ
    analyses), otherwise ENA's FASTA download (e.g. for GCA_937896655)."""
    if record.fasta_ftp:
        ftp = record.fasta_ftp.split(";")[0].strip()
        return ftp if ftp.startswith("http") else "https://" + ftp
    return ENA_FASTA_URL % record.genome_accession


def build_organism_check_candidates(records, studies_by_group_id):
    """Genomes in source groups with 2 or more studies that have no gemsparcl cluster
    yet, each with a download link, so they can be run through gemsparcl query."""
    candidates = []
    for record in sorted(records, key=lambda r: (r.source_group_id, r.genome_accession)):
        if not _multi_study(studies_by_group_id, record.source_group_id):
            continue
        if record.gemsparcl_cluster:
            continue
        candidates.append((record.source_group_id, versioned(record),
                           _fasta_url(record)))
    return candidates


CANDIDATES_FILE = "organism_check_candidates.tsv"
CANDIDATE_COLUMNS = ["source_group_id", "genome_accession", "fasta_url"]


# ---------------------------------------------------------------------------
# Finding duplicates (phase 2)
# ---------------------------------------------------------------------------
#
# Genomes are duplicates when they have the same source samples (same source group), are
# the same organism (same gemsparcl cluster) and come from 2 or more studies.

DUPLICATES_COLUMNS = [
    "duplicate_id", "source_group_id", "gemsparcl_cluster", "n_genomes", "n_studies",
    "declared_runs", "genomes", "studies", "center_names",
]


@dataclass
class DuplicateGroup:
    """All genomes of one organism from one source group, deposited by 2+ studies."""

    duplicate_id: str
    source_group_id: str
    gemsparcl_cluster: str
    genomes: List[str]
    studies: List[str]
    center_names: List[str]
    declared_runs: str


def _compare_declared_runs(genomes):
    """Compare the runs the genomes list in ENA (their declared_runs):

        same       every genome lists runs, and the lists are the same
        different  every genome lists runs, but the lists are not all the same
        unknown    at least one genome lists no runs

    This is for the reader only. It does not decide whether genomes are duplicates,
    because different runs of the same source sample count as the same source data."""
    run_sets = {frozenset(genome.declared_runs) for genome in genomes}
    if frozenset() in run_sets:
        return "unknown"
    return "same" if len(run_sets) == 1 else "different"


def _from_other_study(a, b):
    return bool(a.study_accession and b.study_accession
                and a.study_accession != b.study_accession)


def _organism_verdict(record, members):
    """The duplicate_status of one genome, given every genome in its source group."""
    if not record.gemsparcl_cluster:
        return DUP_NO_CLUSTER
    others = [m for m in members if _from_other_study(record, m)]
    if any(m.gemsparcl_cluster == record.gemsparcl_cluster for m in others):
        return DUP_DUPLICATE
    # A genome from another study with no cluster yet could be the same organism.
    if any(not m.gemsparcl_cluster for m in others):
        return DUP_OTHER_STUDY_NO_CLUSTER
    return DUP_NOT_DUPLICATE


def assign_duplicates(records, studies_by_group_id):
    """Set duplicate_status for every genome in a multi-study group; return the
    duplicate groups."""
    members_by_group = {}
    for record in records:
        if _multi_study(studies_by_group_id, record.source_group_id):
            members_by_group.setdefault(record.source_group_id, []).append(record)

    duplicates = []
    for group_id, members in sorted(members_by_group.items()):
        for record in members:
            record.duplicate_status = _organism_verdict(record, members)

        by_cluster = {}
        for record in members:
            if record.gemsparcl_cluster and record.study_accession:
                by_cluster.setdefault(record.gemsparcl_cluster, []).append(record)
        for gemsparcl_cluster, same in sorted(by_cluster.items()):
            studies = sorted({r.study_accession for r in same})
            if len(studies) < 2:
                continue
            duplicates.append(DuplicateGroup(
                # From the unversioned accessions, so a new version keeps the same id.
                duplicate_id=_short_id("dup_", sorted(r.genome_accession for r in same)),
                source_group_id=group_id,
                gemsparcl_cluster=gemsparcl_cluster,
                genomes=sorted(versioned(r) for r in same),
                studies=studies,
                center_names=sorted({r.center_name for r in same if r.center_name}),
                declared_runs=_compare_declared_runs(same),
            ))
    return duplicates


def write_duplicates_tsv(path, duplicates):
    _write_tsv(path, DUPLICATES_COLUMNS, [
        [dup.duplicate_id, dup.source_group_id, dup.gemsparcl_cluster,
         str(len(dup.genomes)), str(len(dup.studies)), dup.declared_runs,
         _join(dup.genomes), _join(dup.studies), _join(dup.center_names)]
        for dup in sorted(duplicates, key=lambda d: d.duplicate_id)])


# ---------------------------------------------------------------------------
# Command line and running the phases
# ---------------------------------------------------------------------------

def read_accessions(path):
    """Read one accession per line, in order, skipping repeats, blank lines, '#' comments
    and any columns after the first."""
    out, seen = [], set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.split("#", 1)[0]
            token = line.split("\t")[0].split(",")[0].strip()
            if token and token not in seen:
                seen.add(token)
                out.append(token)
    return out


def _apply_trace(record, trace_result):
    """Copy a TraceResult onto a record."""
    if trace_result is not None:
        record.source_samples = list(trace_result.source_samples)
        record.source_sample_status = trace_result.status


def run_pipeline(client, accessions, expand_enabled=True, fasta_dir=None):
    """Phase 1. Look up the input genomes, compare their sequences (with fasta_dir),
    trace them to their source samples, keep one copy of each set of identical genomes,
    find other genomes from the source samples (reverse lookup), trace those too, and
    build the source groups. Returns the records, {source_group_id: studies}, and the
    sets of identical genomes (None without fasta_dir).
    """
    records = lookup_genomes(client, accessions)
    inputs = [records[a.raw] for a in accessions]
    # Every input MAG is checked for identical sequences, including those that can't be
    # traced because ENA lists no sample for them.
    input_mags = [r for r in inputs if not r.source_sample_status]
    identical = find_identical_sequences(input_mags, fasta_dir) if fasta_dir else None
    for record in input_mags:
        if not record.mag_sample:
            record.source_sample_status = STATUS_MAG_SAMPLE_NOT_FOUND
    mags = [r for r in input_mags if not r.source_sample_status]

    sources = trace_source_samples(client, [r.mag_sample for r in mags]) if mags else {}
    for record in mags:
        _apply_trace(record, sources.get(record.mag_sample))
    if identical:
        choose_kept_copies(identical)

    all_records = list(inputs)
    if expand_enabled and sources:
        known = {r.mag_sample for r in inputs if r.mag_sample}
        found = find_other_genomes(client, source_search_keys(sources.values()), known)
        to_trace = [r for r in found if not r.source_sample_status]
        extra_sources = (trace_source_samples(client, [r.mag_sample for r in to_trace])
                         if to_trace else {})
        for record in to_trace:
            _apply_trace(record, extra_sources.get(record.mag_sample))
        all_records.extend(found)

    return all_records, build_source_groups(all_records), identical


def _phase1(args, client):
    LOG.info("")
    LOG.info("Phase 1: provenance of the genomes in ENA")
    accessions = [classify_accession(a) for a in read_accessions(args.accessions)]
    if not accessions:
        LOG.error("no accessions found in %s", args.accessions)
        return 1
    if args.fasta_dir and not os.path.isdir(args.fasta_dir):
        LOG.error("%s is not a folder", args.fasta_dir)
        return 1

    records, studies_by_group_id, identical = run_pipeline(
        client, accessions, expand_enabled=not args.no_expand, fasta_dir=args.fasta_dir)

    write_mags_tsv(os.path.join(args.outdir, "mags.tsv"), records)
    if identical is not None:
        write_identical_sequences_tsv(os.path.join(args.outdir, IDENTICAL_FILE),
                                      identical)

    log_phase1_summary(records, studies_by_group_id, client, identical)
    return 0


def _phase2(args):
    LOG.info("")
    LOG.info("Phase 2: compare gemsparcl clusters within source groups used by 2 or more "
             "studies")
    mags_path = os.path.join(args.outdir, "mags.tsv")
    if not os.path.exists(mags_path):
        LOG.error("%s not found. Run phase 1 first (--accessions).", mags_path)
        return 1

    records = read_mags_tsv(mags_path)
    studies_by_group_id = studies_by_group(records)
    genomes = {versioned(r) for r in records if r.genome_accession}
    attach_clusters(records, read_gemsparcl_clusters(args.gemsparcl_clusters, genomes))

    duplicates = assign_duplicates(records, studies_by_group_id)
    write_duplicates_tsv(os.path.join(args.outdir, "duplicates.tsv"), duplicates)
    write_mags_tsv(mags_path, records)
    candidates = build_organism_check_candidates(records, studies_by_group_id)
    _write_tsv(os.path.join(args.outdir, CANDIDATES_FILE), CANDIDATE_COLUMNS, candidates)

    compared = [r for r in records
                if _multi_study(studies_by_group_id, r.source_group_id)]
    _log_status_table("duplicate_status", compared,
                      (DUP_DUPLICATE, DUP_NOT_DUPLICATE, DUP_OTHER_STUDY_NO_CLUSTER,
                       DUP_NO_CLUSTER))
    _log_row("duplicate groups", len(duplicates),
             "%s genomes, in duplicates.tsv"
             % _number(sum(len(dup.genomes) for dup in duplicates)))
    _log_row("genomes for gemsparcl query", len(candidates),
             "the no_cluster ones, in %s" % CANDIDATES_FILE)
    return 0


HELP_DESCRIPTION = """\
Find MAGs in ENA that were deposited more than once. Two genomes are duplicates when
they come from the same source sample(s), are the same organism (same gemsparcl
cluster) and were deposited by different studies.

phases (give one option to run one phase, or both to run both in one go):
  1. --accessions          look each genome up in ENA, compare their sequences (with
                           --fasta-dir), follow each sample's derived_from links
                           back to the source sample(s), then find every genome built
                           from those samples (reverse lookup)
  2. --gemsparcl-clusters   read the gemsparcl clusters file and report the duplicates"""

HELP_EPILOG = """\
examples:
  # phase 1, then phase 2
  python3 mag_source_duplicates.py --accessions accessions.txt --outdir results/
  python3 mag_source_duplicates.py --gemsparcl-clusters gemsparcl_out/run_clusters.csv \\
      --outdir results/

  # both phases in one go
  python3 mag_source_duplicates.py --accessions accessions.txt \\
      --gemsparcl-clusters gemsparcl_out/run_clusters.csv --outdir results/

output files, in --outdir:
  mags.tsv                       one row per genome: source samples, study,
                                 source_sample_status, duplicate_status, sequence_check
  identical_sequences.tsv        with --fasta-dir: one row per set of input genomes
                                 with identical sequences
  organism_check_candidates.tsv  genomes with no gemsparcl cluster, with a FASTA link;
                                 after phase 2 they can be run through gemsparcl query
  duplicates.tsv                 phase 2: one row per duplicate group
  mag_source_duplicates.log      what was printed on screen, for every run"""


def _help_formatter(prog):
    """Keep the description and examples as written, give option names room so their
    help starts on the same line, and wrap at 90 characters whatever the terminal width."""
    return argparse.RawDescriptionHelpFormatter(prog, max_help_position=30, width=90)


def build_arg_parser():
    parser = argparse.ArgumentParser(
        prog="mag_source_duplicates.py", description=HELP_DESCRIPTION, epilog=HELP_EPILOG,
        formatter_class=_help_formatter)
    parser.add_argument("--accessions", metavar="FILE",
                        help="run phase 1 on this file of genome accessions (GCA_, ERZ or "
                             "WGS set), one per line")
    parser.add_argument("--gemsparcl-clusters", metavar="FILE",
                        help="run phase 2 with this gemsparcl clusters file "
                             "(*_clusters.csv, columns genome_id and cluster)")
    parser.add_argument("--outdir", metavar="DIR", default="results",
                        help="folder for the output files; phase 2 reads mags.tsv from "
                             "here (default: %(default)s)")
    parser.add_argument("--fasta-dir", metavar="DIR",
                        help="phase 1: compare the input genomes' FASTA files in this "
                             "folder (named by accession, e.g. GCA_937896655.1.fna.gz); "
                             "of genomes with identical sequences, only one is "
                             "compared further")
    parser.add_argument("--cache-dir", metavar="DIR",
                        help="keep ENA answers in this folder, so a rerun reuses them "
                             "instead of asking ENA again")
    parser.add_argument("--no-expand", action="store_true",
                        help="skip the reverse lookup (finding other genomes from the "
                             "same source samples); only the input genomes are compared "
                             "(finds far fewer duplicates)")
    return parser


class _MessageFormatter(logging.Formatter):
    """Show only the message; warnings and errors also start with their level."""

    def format(self, record):
        message = record.getMessage()
        if record.levelno >= logging.WARNING:
            return "%s: %s" % (record.levelname, message)
        return message


def _add_log_handler(handler, level):
    handler.setLevel(level)
    handler.setFormatter(_MessageFormatter())
    LOG.addHandler(handler)


def _start_logging():
    """Print the summary, warnings and errors on screen. The log file is added once
    --outdir exists."""
    for handler in list(LOG.handlers):  # so that calling main() twice doesn't print twice
        LOG.removeHandler(handler)
        handler.close()
    LOG.setLevel(logging.DEBUG)
    LOG.propagate = False
    _add_log_handler(logging.StreamHandler(sys.stderr), logging.INFO)


def _start_log_file(path):
    """Add every message, the step-by-step ones too, to the log file. Each run is
    appended, because phase 1 and phase 2 are often run separately into one folder."""
    if os.path.exists(path) and os.path.getsize(path):
        with open(path, "a", encoding="utf-8") as handle:
            handle.write("\n")  # a blank line between runs
    _add_log_handler(logging.FileHandler(path, encoding="utf-8"), logging.DEBUG)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    args = build_arg_parser().parse_args(argv)
    _start_logging()
    if not args.accessions and not args.gemsparcl_clusters:
        LOG.error("give --accessions (phase 1) or --gemsparcl-clusters (phase 2)")
        return 2
    try:
        os.makedirs(args.outdir, exist_ok=True)
        log_path = os.path.join(args.outdir, LOG_FILE)
        _start_log_file(log_path)
        LOG.info("mag_source_duplicates.py, %s", time.strftime("%Y-%m-%d %H:%M"))
        LOG.info("command: %s", " ".join(shlex.quote(a) for a in [sys.argv[0]] + argv))
        code = 0
        if args.accessions:
            client = EnaClient(cache_dir=args.cache_dir)
            code = _phase1(args, client)
        if code == 0 and args.gemsparcl_clusters:
            code = _phase2(args)
        if code == 0:
            LOG.info("")
            LOG.info("Results and %s written to %s", LOG_FILE, args.outdir)
        return code
    except (EnaError, ValueError, OSError) as exc:
        LOG.error("%s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
