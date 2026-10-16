#!/usr/bin/env python3
"""
Identify MAGs independently generated from the same source data and deposited in ENA.

Input:  a file of INSDC genome accessions (GCA_*, ERZ*, or a WGS set e.g. CABKLP010000000).
Output: mags.tsv (one row per genome), source_groups.tsv (one row per source group,
i.e. one set of raw reads), organism_check_candidates.tsv (genomes whose organism is
still unknown), and duplicates.tsv after the organism join. A source group is NOT a
gemsparcl cluster/GCU.

A duplicate is the conjunction of three claims, and no two of them suffice:

    1. same source data      the genomes trace to the same source sample(s)
    2. different submitters   they come from different ENA studies
    3. same organism          same GCU cluster reported from gemsparcl
"""

import argparse
import csv
import hashlib
import io
import json
import logging
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set

LOG = logging.getLogger("find_mag_source_duplicates")

# ---------------------------------------------------------------------------
# Input genome accession kinds
# ---------------------------------------------------------------------------

KIND_GCA = "GCA"
KIND_ERZ = "ERZ"
KIND_WGS = "WGS"
KIND_INVALID = "invalid"

# Accession patterns, taken from the Portal API's own rejection message:
#     wgs  ^[A-Z]{4}[A-Z]{2}?[0-9]{2}[0]{6,9}$
# A WGS set accession has a 4- or 6-letter prefix, two digits of version, then zeros. ENA
# rejects the unpadded short form, so we pad it.
RE_GCA = re.compile(r"^GC[AF]_[0-9]{9}$")
RE_ERZ = re.compile(r"^ERZ[0-9]+$")
RE_WGS_PADDED = re.compile(r"^[A-Z]{4}([A-Z]{2})?[0-9]{2}0{6,9}$")
RE_WGS_SHORT = re.compile(r"^[A-Z]{4}([A-Z]{2})?[0-9]{2}$")

# A real version suffix: a dot, digits, end of string. Nothing else counts.
RE_VERSION_SUFFIX = re.compile(r"^(.*)\.[0-9]+$")

WGS_PAD = "0" * 7


@dataclass(frozen=True)
class Accession:
    """
    One input accession: what the user wrote, and the form we query ENA with.
    `raw` is kept so every output row traces back to exactly what was supplied, even when
    the queried form differs (a stripped version suffix, or a padded WGS set).
    """

    raw: str
    kind: str
    query: str


def classify_accession(raw):
    """Classify and normalise one input accession.

    The version suffix is stripped, because ENA ignores it and returns every version as a
    row; we pick the latest ourselves. Short WGS accessions are zero-padded, as ENA rejects
    `CABKLP01` but accepts `CABKLP010000000`. Anything unrecognised returns KIND_INVALID
    rather than raising, so one bad line cannot abort a run.
    """
    original = raw.strip() if raw else ""
    token = original.upper()
    if not token:
        return Accession(raw=original, kind=KIND_INVALID, query="")

    version_match = RE_VERSION_SUFFIX.match(token)
    unversioned = version_match.group(1) if version_match else token

    if RE_GCA.match(unversioned):
        return Accession(raw=original, kind=KIND_GCA, query=unversioned)
    if RE_ERZ.match(unversioned):
        return Accession(raw=original, kind=KIND_ERZ, query=unversioned)
    if RE_WGS_PADDED.match(unversioned):
        return Accession(raw=original, kind=KIND_WGS, query=unversioned)
    if RE_WGS_SHORT.match(unversioned):
        return Accession(raw=original, kind=KIND_WGS, query=unversioned + WGS_PAD)
    return Accession(raw=original, kind=KIND_INVALID, query="")


# ---------------------------------------------------------------------------
# The "I was made from X" field
# ---------------------------------------------------------------------------
#
# `related_sample_accession` on the MAG's own sample record. Submitter-typed free text, so
# four shapes appear, and mishandling any of them loses duplicates silently:
#
#   biosample   SAMEA2622837           the ordinary case
#   secondary   ERS490230              the SAME sample, other spelling. ENA does not
#                                      cross-index them, so one form misses whole cohorts
#   run         ERR982800              a RUN. 16-100% of values, submitter-dependent.
#                                      Resolve to its sample or cohorts split apart
#   range       ERS9142775-ERS9142786  a block of samples. NOT expanded: the submitter's
#                                      own explicit list omits a member it covers
#
# Also ';'-separated and multi-valued (co-assembly), in no fixed order.

PARENT_BIOSAMPLE = "biosample"
PARENT_SECONDARY = "secondary"
PARENT_RUN = "run"
PARENT_RANGE = "range"
PARENT_UNPARSABLE = "unparsable"

RELATION_DERIVED_FROM = "derived_from"
RELATION_COMPOSED_OF = "composed_of"

RE_BIOSAMPLE = re.compile(r"^SAM[END]A?[0-9]{6,}$")
RE_SECONDARY = re.compile(r"^[EDSR]RS[0-9]{6,}$")
RE_RUN = re.compile(r"^[EDSR]RR[0-9]{6,}$")
RE_RANGE = re.compile(r"^[EDSR]RS[0-9]{6,}-[EDSR]RS[0-9]{6,}$")

# For scraping the run of an ERZ bin out of free text, which is the only place it is
# recorded. Nothing else is scraped: an accession appearing in prose proves co-occurrence,
# never a relationship, so prose cannot establish provenance.
RE_RUN_IN_TEXT = re.compile(r"[EDSR]RR[0-9]{6,}")

_STRIP_CHARS = ' \t"\''


@dataclass(frozen=True)
class ParentValue:
    """One parsed `X:relation` entry: the accession, and what kind of thing it is."""

    raw: str
    kind: str


def classify_parent(token):
    """Work out which of the four known shapes a parent value is."""

    if RE_BIOSAMPLE.match(token):
        return PARENT_BIOSAMPLE
    if RE_SECONDARY.match(token):
        return PARENT_SECONDARY
    if RE_RUN.match(token):
        return PARENT_RUN
    if RE_RANGE.match(token):
        return PARENT_RANGE
    return PARENT_UNPARSABLE


def parse_related_field(value, relation=RELATION_DERIVED_FROM):
    """Parse `related_sample_accession` into the parents carrying `relation`.

    Returns EVERY matching value, in the order ENA stored them: a co-assembled MAG lists
    several parents, and sorting is left to the fingerprint so each row's parent kinds stay
    lined up. Unrecognised shapes are returned flagged, never dropped.
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
        if found_relation.strip(_STRIP_CHARS) != relation or not token:
            continue
        out.append(ParentValue(raw=token, kind=classify_parent(token)))
    return out


def _dedupe(items):
    """Order-preserving de-duplication, so output stays deterministic."""
    seen = set()
    out = []
    for item in items:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out



def extract_runs_from_text(text):
    """Scrape RUN accessions out of free text.

    For ERZ-level bins the structured field is always empty, so this is the only route to a
    run. Evidence only, never a grouping key.
    """
    if not text:
        return []
    return _dedupe(RE_RUN_IN_TEXT.findall(text))



# ---------------------------------------------------------------------------
# Talking to ENA
# ---------------------------------------------------------------------------
#
# Four requirements:
#
#   POST, not GET   a GET with ~120 OR-terms is rejected as too long; large GET accession
#                   batches time out
#   Check the BODY  a transient failure arrives as HTTP 200 with "ERROR occurred" in the
#                   body. Trusting the status code turns a glitch into "no duplicates"
#   Split by form   one request may not mix SAM* and ERS*/SRS* accessions; it 400s
#   Cache to disk   an HPC job killed at the walltime boundary restarts cheaply
#
# Hence one rule: a failure RAISES, never returns an empty list. "No results" and "the
# request failed" sharing a code path produced 81 false findings in our own checker.
#
# Standard library only -- nothing to install on an HPC node.

ENA_SEARCH_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_ERROR_MARKER = "ERROR occurred"

# HTTP codes that mean "this request is wrong and always will be", so retrying is waste:
# a 400 is what ENA returns for an accession sent to the wrong result type or for
# mixed accession forms.
PERMANENT_HTTP_CODES = frozenset({400, 404, 414})


class EnaError(RuntimeError):
    """A request could not be completed. Never raised to mean "nothing matched"."""


class EnaHttpError(Exception):
    """An HTTP error from ENA, carrying the status code so we can judge retryability."""

    def __init__(self, status_code, body=""):
        super().__init__("HTTP %s" % status_code)
        self.status_code = status_code
        self.body = body


class _Response:
    """The little bit of a response this client needs."""

    def __init__(self, text, status_code=200):
        self.text = text
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise EnaHttpError(self.status_code, self.text)


class _UrllibSession:
    """Default transport: POST over the standard library. No third-party dependency."""

    def post(self, url, data=None, timeout=None):
        encoded = urllib.parse.urlencode(data or {}).encode("utf-8")
        request = urllib.request.Request(url, data=encoded)
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return _Response(response.read().decode("utf-8", errors="replace"))
        except urllib.error.HTTPError as exc:
            body = ""
            try:
                body = exc.read().decode("utf-8", errors="replace")
            except Exception:  # pragma: no cover - body already consumed
                pass
            raise EnaHttpError(exc.code, body)


def split_accessions_by_form(accessions):
    """Group accessions so that no request mixes the two sample-accession forms.

    ENA rejects a batch containing both `SAMEA…` and `ERS…` style accessions with HTTP
    400. Anything that is neither form (genome accessions, runs) gets its own group.
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
    """Fetches records from the ENA Portal API."""

    def __init__(self, session=None, cache_dir=None, timeout=180, max_retries=4,
                 retry_wait=5.0, batch_size=500):
        self.session = session if session is not None else _UrllibSession()
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
        if not self.cache_dir:
            return None
        key = json.dumps(payload, sort_keys=True)
        digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
        return os.path.join(self.cache_dir, digest + ".tsv")

    def _post(self, payload):
        """One request, with retries. Returns the body text, or raises EnaError."""
        path = self._cache_path(payload)
        if path and os.path.exists(path):
            self.cache_hits += 1
            with open(path, "r", encoding="utf-8") as handle:
                return handle.read()

        last_error = None
        for attempt in range(1, self.max_retries + 1):
            try:
                self.request_count += 1
                response = self.session.post(ENA_SEARCH_URL, data=payload,
                                             timeout=self.timeout)
                response.raise_for_status()
                text = response.text
                if ENA_ERROR_MARKER in text:
                    # HTTP 200 with the error inside the body.
                    last_error = "ENA returned its in-body error marker"
                    LOG.warning("ENA in-body error, attempt %s/%s (result=%s)",
                                attempt, self.max_retries, payload.get("result"))
                else:
                    if path:
                        with open(path, "w", encoding="utf-8") as handle:
                            handle.write(text)
                    return text
            except EnaHttpError as exc:
                last_error = "HTTP %s %s" % (exc.status_code, exc.body[:200])
                if exc.status_code in PERMANENT_HTTP_CODES:
                    # Wrong request. Retrying cannot help.
                    raise EnaError(
                        "ENA rejected the request (result=%s): %s"
                        % (payload.get("result"), last_error)
                    )
                LOG.warning("ENA HTTP %s, attempt %s/%s", exc.status_code, attempt,
                            self.max_retries)
            except Exception as exc:
                last_error = str(exc)
                LOG.warning("ENA request failed, attempt %s/%s: %s", attempt,
                            self.max_retries, exc)
            if attempt < self.max_retries and self.retry_wait:
                time.sleep(self.retry_wait * attempt)

        raise EnaError("ENA request failed after %s attempts (result=%s): %s"
                       % (self.max_retries, payload.get("result"), last_error))

    @staticmethod
    def _parse_tsv(text):
        """Turn a TSV body into dictionaries, padding short rows rather than losing them."""
        if not text.strip():
            return []
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

    def search(self, result, fields, query=None, accessions=None, limit=0):
        """Run one search and return the rows.

        Give either `query` or `accessions`. Accession lists are de-duplicated, grouped so
        no request mixes accession forms, and chunked. `limit=0` means everything and
        is not silently capped.
        """
        base = {
            "result": result,
            "fields": ",".join(fields),
            "format": "tsv",
            "limit": str(limit),
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
# Resolving an accession to a genome record
# ---------------------------------------------------------------------------

STATUS_RESOLVED = "resolved"
STATUS_SAMPLE_ONLY = "sample_only"
STATUS_NO_DERIVED_FROM = "unresolved_no_derived_from"
STATUS_NOT_FOUND = "unresolved_not_found"
STATUS_INVALID = "invalid_accession"
STATUS_NOT_A_MAG = "not_a_mag"
STATUS_UNPARSABLE_PARENT = "unresolved_unparsable_parent"
STATUS_PARENT_NOT_PUBLIC = "unresolved_parent_not_public"
# Deliberately invalid: the walk has four exit branches, and a real default would make a
# missed one report "the submitter recorded nothing" instead of a bug.
STATUS_UNSET = "BUG_status_never_set"

NON_MAG_ASSEMBLY_TYPES = frozenset({
    "clone or isolate",                                 # cultured; no metagenome source
    "environmental single-cell amplified genome (sag)",  # one cell, not assembled from a
                                                        # metagenome, so nothing to share
    "primary metagenome",                               # the whole assembly, not a genome
})

# A SEPARATE field list per result type. One shared list is rejected outright -- the three
# endpoints offer genuinely different fields, and the field index is neither exhaustive nor
# reliable, so these were verified by querying.
#
# Absences that matter:
#   `broker_name`        absent on assembly/wgs_set, and restates center_name
#   `first_public`       absent on assembly, so last_updated is the only date for the
#                        2.5M GCA-level MAGs
#   `sample_description` absent on assembly/wgs_set; `description` is the universal one
#   wgs_set              has no centre, software or quality fields at all -- which is why
#                        a WGS input is resolved via its GCA instead

ASSEMBLY_FIELDS = [
    "accession", "version", "study_accession", "center_name", "sample_accession",
    "run_accession", "binning_software", "assembly_software", "completeness_score",
    "contamination_score", "assembly_type", "last_updated", "description", "sample_title",
    "tax_id", "scientific_name",
]
ANALYSIS_FIELDS = [
    "analysis_accession", "first_public", "study_accession", "center_name",
    "sample_accession", "run_accession", "binning_software", "assembly_software",
    "completeness_score", "contamination_score", "assembly_type", "last_updated",
    "description", "sample_title", "sample_description", "tax_id", "scientific_name",
    # ERZ analyses are not served by the browser FASTA endpoint; this is where
    # their contig.fa.gz lives.
    "generated_ftp",
]
WGS_FIELDS = [
    "accession", "assembly_accession", "first_public", "study_accession",
    "sample_accession", "run_accession", "assembly_type", "last_updated", "description",
    "tax_id", "scientific_name",
]

RUN_EVIDENCE_STRUCTURED = "structured"
RUN_EVIDENCE_DESCRIPTION = "description"
RUN_EVIDENCE_NONE = "none"


@dataclass
class GenomeRecord:
    """One genome, as resolved from ENA, plus what later stages work out about it."""

    accession: str
    input_accession: str
    genome_type: str
    resolution_status: str
    version: str = ""
    mag_sample: str = ""
    study_accession: str = ""
    center_name: str = ""
    last_updated: str = ""
    first_public: str = ""
    assembly_software: str = ""
    binning_software: str = ""
    completeness_score: str = ""
    contamination_score: str = ""
    assembly_type: str = ""
    assembly_accession: str = ""
    tax_id: str = ""
    scientific_name: str = ""
    generated_ftp: str = ""
    wgs_set: str = ""
    description_text: str = ""
    declared_runs: List[str] = field(default_factory=list)
    run_evidence: str = RUN_EVIDENCE_NONE
    in_input: bool = True
    # Filled in by later stages.
    source_samples: List[str] = field(default_factory=list)
    source_parent_types: List[str] = field(default_factory=list)
    # The runs of the SOURCE SAMPLE, from the walk. A property of the group, identical for
    # every member, so it must NOT be used to compare members -- see declared_runs.
    source_sample_runs: List[str] = field(default_factory=list)
    # True when SOME parents resolved and others did not, so the source set is a known
    # SUBSET of the truth: kept, flagged, and barred from creating a generation event.
    source_partial: bool = False
    source_group_id: str = ""
    relationship: str = ""
    # Organism identity, from gemsparcl.
    gcu_id: str = ""
    gcu_note: str = ""


# Which result type answers for which kind of accession. There is exactly one correct
# endpoint per kind: sending an accession to the wrong one returns HTTP 400, not an empty
# answer, so a try-each-in-turn cascade cannot work.
_ROUTES = (
    (KIND_GCA, "assembly", ASSEMBLY_FIELDS, "accession"),
    (KIND_ERZ, "analysis", ANALYSIS_FIELDS, "analysis_accession"),
    (KIND_WGS, "wgs_set", WGS_FIELDS, "accession"),
)


# GenomeRecord fields whose name is also the ENA field name, copied straight across.
_PASSTHROUGH_FIELDS = (
    "version", "study_accession", "center_name", "last_updated", "first_public",
    "assembly_software", "binning_software", "completeness_score",
    "contamination_score", "assembly_type", "assembly_accession",
    "tax_id", "scientific_name", "generated_ftp",
)


def _row_to_record(row, acc_key, input_accession, genome_type):
    """Turn one ENA row into a GenomeRecord, working out the run evidence as we go."""
    description_text = " ".join([row.get("description", ""),
                                 row.get("sample_description", ""),
                                 row.get("sample_title", "")]).strip()

    declared = [r for r in re.split(r"[;,]", row.get("run_accession", "")) if r.strip()]
    if declared:
        evidence = RUN_EVIDENCE_STRUCTURED
    else:
        # ERZ-level bins never carry a structured run, so the description is the only
        # route. Evidence only -- never used as a grouping key.
        declared = extract_runs_from_text(description_text)
        evidence = RUN_EVIDENCE_DESCRIPTION if declared else RUN_EVIDENCE_NONE

    return GenomeRecord(
        accession=row.get(acc_key, ""),
        input_accession=input_accession,
        genome_type=genome_type,
        resolution_status=STATUS_RESOLVED,
        mag_sample=row.get("sample_accession", ""),
        description_text=description_text,
        declared_runs=[r.strip() for r in declared],
        run_evidence=evidence,
        **{name: row.get(name, "") for name in _PASSTHROUGH_FIELDS}
    )


def _version_key(value):
    """Versions sort numerically: GCA_000001405 reached version 29, and "29" < "9"."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return -1


def _latest(items, get_version):
    """Several rows for one accession means several versions."""
    return max(items, key=lambda item: _version_key(get_version(item)))


def resolve_genomes(client, accessions):
    """Resolve input accessions to genome records, keyed by what the user typed.

    Every input gets a record. An accession that cannot be resolved comes back with a
    `resolution_status` saying why, rather than being dropped, so nothing disappears
    between the input file and the output.
    """
    records = {}
    wanted = []
    for accession in accessions:
        if accession.kind == KIND_INVALID:
            records[accession.raw] = GenomeRecord(
                accession=accession.raw,
                input_accession=accession.raw,
                genome_type=KIND_INVALID,
                resolution_status=STATUS_INVALID,
            )
        else:
            wanted.append(accession)

    if not wanted:
        return records

    by_query = {a.query: a for a in wanted}

    for kind, result, fields, acc_key in _ROUTES:
        queries = sorted(a.query for a in wanted if a.kind == kind)
        if not queries:
            continue
        candidates = {}
        for row in client.search(result, fields, accessions=queries):
            key = row.get(acc_key, "")
            if key in by_query:
                candidates.setdefault(key, []).append(
                    _row_to_record(row, acc_key, by_query[key].raw, kind))
        for key, found in candidates.items():
            records[by_query[key].raw] = _latest(found, lambda r: r.version)

    # A WGS set carries almost no evidence of its own -- no centre, no software, no quality
    # scores. It does carry its GCA, so re-resolve through `assembly` to
    # get a full record. Without this, WGS-supplied rows would have blank evidence columns
    # and their independence call could not be audited.
    wgs_to_gca = {}
    for accession in wanted:
        if accession.kind != KIND_WGS:
            continue
        record = records.get(accession.raw)
        if record is None or record.resolution_status != STATUS_RESOLVED:
            continue
        gca = (record.assembly_accession or "").strip()
        if gca:
            wgs_to_gca[accession.raw] = gca
    if wgs_to_gca:
        upgraded = {}
        for row in client.search("assembly", ASSEMBLY_FIELDS,
                                 accessions=sorted(set(wgs_to_gca.values()))):
            key = row.get("accession", "")
            if key:
                upgraded.setdefault(key, []).append(row)
        for raw, gca in wgs_to_gca.items():
            rows = upgraded.get(gca)
            if not rows:
                continue
            best = _latest(rows, lambda r: r.get("version"))
            full = _row_to_record(best, "accession", raw, KIND_WGS)
            # Keep the WGS accession the user supplied as the identity of the row.
            full.wgs_set = records[raw].accession
            full.assembly_accession = gca
            records[raw] = full

    # Anything requested and not returned simply has no record: a valid-shape accession
    # with nothing behind it returns 200 and zero rows, which is distinguishable from the
    # wrong-endpoint 400.
    for accession in wanted:
        if accession.raw not in records:
            records[accession.raw] = GenomeRecord(
                accession=accession.query,
                input_accession=accession.raw,
                genome_type=accession.kind,
                resolution_status=STATUS_NOT_FOUND,
            )

    for record in records.values():
        if record.resolution_status != STATUS_RESOLVED:
            continue
        if (record.assembly_type or "").lower() in NON_MAG_ASSEMBLY_TYPES:
            record.resolution_status = STATUS_NOT_A_MAG

    return records


# ---------------------------------------------------------------------------
# Walking up to the source data
# ---------------------------------------------------------------------------
#
# A MAG has its own registered, virtual sample, which points back at the real material
# with `derived_from`. Following that pointer is this whole stage, and every way of getting
# it wrong is silent: the run completes and reports a differently-shaped answer.
#
# The walk ends when a sample HAS READS. One `read_run` lookup answers three questions at
# once -- has it reads, which ones, and what is its other accession form -- so that is the
# terminating test. Batched per depth level.

SOURCE_PARENT_SAMPLE = "sample"
SOURCE_PARENT_RUN_RESOLVED = "run_resolved"
SOURCE_PARENT_RANGE_LITERAL = "range_literal"

SAMPLE_FIELDS = ["sample_accession", "secondary_sample_accession",
                 "related_sample_accession", "description", "sample_title"]
READ_RUN_FIELDS = ["run_accession", "sample_accession", "secondary_sample_accession"]

# Trails of two steps occur; ten is a ceiling, not an expectation. The per-MAG visited set
# already stops cycles, but a malformed loop must not cost a whole HPC job.
MAX_WALK_DEPTH = 10


@dataclass
class SourceResolution:
    """Where one MAG sample's source data was traced to."""

    source_samples: List[str] = field(default_factory=list)   # canonical, order found
    parent_types: List[str] = field(default_factory=list)     # aligned with source_samples
    forms: Set[str] = field(default_factory=set)              # every form, for expansion
    runs: List[str] = field(default_factory=list)
    status: str = STATUS_UNSET
    # Some parents resolved, others did not. The source set is then a known subset of the
    # truth, so the group key is wrong and the MAG groups alone, away from its real
    # siblings. Kept in the output, flagged, barred from the independence call.
    partial: bool = False


def _both_forms(row):
    """The two accessions ENA gives for one sample: BioSample, and INSDC secondary."""
    return [row.get("sample_accession", "").strip(),
            row.get("secondary_sample_accession", "").strip()]


def _canonical_sample(row, fallback=""):
    """The form used as the source-group key: BioSample where known."""
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
    """Index `read_run` rows by run accession and by either form of their sample."""
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


class _Walk:
    """One MAG sample's walk: what it found, what is left to check, where it has been.

    State is per MAG sample even though requests are shared: two MAGs routinely share a
    parent, and a shared visited-set would rob the second one of its source.
    """

    def __init__(self, mag_sample):
        self.result = SourceResolution()
        self.visited = {mag_sample}
        self.frontier = []
        self.saw_unparsable = False
        self.saw_missing_parent = False

    def add_source(self, canonical, parent_type, forms=(), runs=()):
        if canonical and canonical not in self.result.source_samples:
            self.result.source_samples.append(canonical)
            self.result.parent_types.append(parent_type)
        self.result.forms.update(form for form in forms if form)
        for run in runs:
            if run and run not in self.result.runs:
                self.result.runs.append(run)

    def enqueue(self, parent):
        """Act on one parsed parent value according to its shape."""
        if parent.kind == PARENT_RANGE:
            # Never expanded, so the range text itself is the identity. These MAGs group
            # with each other and nothing else, because we cannot prove a range and an
            # explicit list mean the same samples.
            self.add_source(parent.raw, SOURCE_PARENT_RANGE_LITERAL, forms=[parent.raw])
        elif parent.kind == PARENT_UNPARSABLE:
            self.saw_unparsable = True
        elif parent.raw not in self.visited:
            self.visited.add(parent.raw)
            self.frontier.append(parent)

    def finish(self):
        # Finding any source outranks the failures. But a failure alongside a success
        # means the set is incomplete, which makes the group key wrong, not approximate.
        if self.result.source_samples:
            self.result.status = STATUS_RESOLVED
            self.result.partial = self.saw_unparsable or self.saw_missing_parent
        elif self.saw_unparsable:
            self.result.status = STATUS_UNPARSABLE_PARENT
        elif self.saw_missing_parent:
            self.result.status = STATUS_PARENT_NOT_PUBLIC
        else:
            self.result.status = STATUS_NO_DERIVED_FROM
        return self.result


def resolve_sources(client, mag_samples):
    """Trace each MAG sample to its source sample(s).

    Every requested sample gets an entry, with a status saying why when nothing was found.
    A MAG whose provenance is unrecoverable still has to appear in the output.
    """
    wanted = _dedupe(s.strip() for s in mag_samples if s and s.strip())
    if not wanted:
        return {}

    walks = {sample: _Walk(sample) for sample in wanted}

    # One request for all the MAG samples, not one each.
    index = _index_samples_by_either_form(
        client.search("sample", SAMPLE_FIELDS, accessions=wanted))
    for mag_sample, walk in walks.items():
        row = index.get(mag_sample)
        if row is None:
            # No record for the MAG's own sample, so there is nothing to walk.
            continue
        walk.visited.update(form for form in _both_forms(row) if form)
        for parent in parse_related_field(row.get("related_sample_accession")):
            walk.enqueue(parent)

    for depth in range(1, MAX_WALK_DEPTH + 1):
        frontier = {}
        for walk in walks.values():
            for parent in walk.frontier:
                frontier[parent.raw] = parent.kind
        if not frontier:
            break
        LOG.info("source walk depth %s: %s candidate parent(s)", depth, len(frontier))

        by_run, by_sample = _index_runs(
            client.search("read_run", READ_RUN_FIELDS, accessions=sorted(frontier)))

        read_less = {}
        for mag_sample, walk in walks.items():
            current, walk.frontier = walk.frontier, []
            for parent in current:
                if parent.kind == PARENT_RUN:
                    row = by_run.get(parent.raw)
                    if row is None or not _canonical_sample(row):
                        # A run cited as a parent that ENA will not serve, and a run has
                        # no parent of its own, so the trail ends here.
                        walk.saw_missing_parent = True
                        continue
                    walk.add_source(_canonical_sample(row), SOURCE_PARENT_RUN_RESOLVED,
                                    forms=_both_forms(row), runs=[parent.raw])
                    continue
                rows = by_sample.get(parent.raw)
                if rows:
                    walk.add_source(
                        _canonical_sample(rows[0], fallback=parent.raw),
                        SOURCE_PARENT_SAMPLE,
                        forms=_both_forms(rows[0]) + [parent.raw],
                        runs=[r.get("run_accession", "").strip() for r in rows],
                    )
                else:
                    read_less.setdefault(parent.raw, []).append(mag_sample)

        if not read_less:
            continue
        # A parent with no reads is either an intermediate bin with a parent of its own,
        # or a reference ENA will not serve. One request decides which.
        parents = _index_samples_by_either_form(
            client.search("sample", SAMPLE_FIELDS, accessions=sorted(read_less)))
        for token in sorted(read_less):
            row = parents.get(token)
            for mag_sample in read_less[token]:
                walk = walks[mag_sample]
                if row is None:
                    walk.saw_missing_parent = True
                    continue
                walk.visited.update(form for form in _both_forms(row) if form)
                for parent in parse_related_field(row.get("related_sample_accession")):
                    walk.enqueue(parent)

    unfinished = sum(1 for walk in walks.values() if walk.frontier)
    if unfinished:
        LOG.warning("source walk hit the %s-step ceiling with %s trail(s) still open",
                    MAX_WALK_DEPTH, unfinished)

    return {mag_sample: walk.finish() for mag_sample, walk in walks.items()}


# ---------------------------------------------------------------------------
# Reverse expansion: finding the depositors who are NOT in the input list
# ---------------------------------------------------------------------------
#
# If a MAG in the input is a re-analysis of someone else's work, the other group's MAGs
# are by definition not in that input. Without this reverse lookup the script could only
# find pairs the user already submitted together.
#
# Four key types per source, or the miss is DIRECTIONAL: you find a cohort only if you
# happened to list one of its MAGs.
#
#   biosample form   SAMEA3504092           }  DISJOINT result sets; in 2 of 7 measured
#   secondary form   ERS811241              }  samples one of them returns ZERO
#   every run        ERR982800                 without it, run-declaring depositors are
#                                              invisible -- 3 centres read as 2
#   range literal    ERS9142775-ERS9142786     a range-declarer is not found by querying
#                                              one of its members

DEFAULT_OR_BATCH = 100


def build_expansion_query(keys):
    """OR-join reverse-lookup terms into one `query=` expression.

    Mixing the two accession forms here is SAFE, unlike in an `includeAccessions` batch,
    which 400s.
    """
    terms = ['related_sample_accession="%s:derived_from"' % key for key in keys if key]
    return " OR ".join(terms)


def expansion_keys(resolutions):
    """Every key worth reverse-querying, deduplicated and sorted.

    A PARTIAL resolution still contributes: the parents it did find are real source data.
    `partial` restricts the independence call, not what we may look for.
    """
    keys = set()
    for resolution in resolutions:
        if resolution.status != STATUS_RESOLVED:
            continue
        keys.update(form for form in resolution.forms if form)
        keys.update(run for run in resolution.runs if run)
    return sorted(keys)


def expand(client, keys, or_batch_size=DEFAULT_OR_BATCH):
    """Reverse-lookup every key; return the discovered MAG sample accessions, in order.

    A failed request RAISES rather than returning nothing: an expansion that silently
    yields [] reads as "no duplicates exist".
    """
    keys = [key for key in keys if key]
    found = []
    seen = set()
    for start in range(0, len(keys), or_batch_size):
        query = build_expansion_query(keys[start:start + or_batch_size])
        if not query:
            continue
        for row in client.search("sample", SAMPLE_FIELDS, query=query):
            # Take whichever accession form ENA returned, so a sample is never dropped for
            # having been keyed under the other one.
            accession = _canonical_sample(row)
            if accession and accession not in seen:
                seen.add(accession)
                found.append(accession)
    LOG.info("expansion: %s key(s) -> %s discovered sample(s)", len(keys), len(found))
    return found


# ---------------------------------------------------------------------------
# Source groups, independence, and the organism-check gate
# ---------------------------------------------------------------------------
#
# Grouping is on the SORTED source-sample set, because both orderings occur within one
# study. Nothing here may say "duplicate": organism identity is unknown until the join.

# Phase 1 verdicts.
REL_PENDING_ORGANISM_CHECK = "pending_organism_check"
REL_SAME_STUDY_REDUNDANCY = "same_study_redundancy"
REL_SINGLETON = "singleton"
REL_UNRESOLVED = "unresolved"
REL_PARTIAL_SOURCE = "partial_source"


@dataclass
class SourceGroup:
    """Genomes that declare the SAME SOURCE SAMPLE(S).

    Same source *sample*, not same reads: so two groups could have used different subsets.
    Reported as `run_mismatch`.

    NOT a gemsparcl cluster and NOT a duplicate set -- the members are usually different
    organisms.
    """

    source_group_id: str
    source_samples: List[str]
    n_genomes: int = 0
    n_events: int = 0
    independent: bool = False
    needs_organism_check: bool = False
    n_distinct_run_sets: int = 0
    run_coverage: float = 0.0
    # Members declared different runs. A FLAG, never a verdict: two runs of one library
    # are technical replicates.
    run_mismatch: bool = False
    n_ineligible: int = 0
    studies: List[str] = field(default_factory=list)
    centres: List[str] = field(default_factory=list)


def fingerprint(source_samples):
    """Stable group id from the SORTED source-sample set.

    ENA stores multi-value parents in no fixed order, and the same two samples appear both
    ways round within one study, so keying on the raw string would split one submission.
    """
    items = sorted(set(sample for sample in source_samples if sample))
    digest = hashlib.sha1("|".join(items).encode("utf-8")).hexdigest()
    return "sg_" + digest[:12]



def _is_eligible(record):
    """May this row count toward the independence call?

    No if its source set is a known subset of the truth: the group key is then wrong. The
    row stays in the output, but cannot create a generation event or trigger a download.
    """
    return not record.source_partial


def build_source_groups(records):
    """Assign group ids, compute stats, set verdicts, and gate the organism check."""
    grouped = {}
    for record in records:
        if record.resolution_status != STATUS_RESOLVED or not record.source_samples:
            record.source_group_id = ""
            record.relationship = REL_UNRESOLVED
            continue
        record.source_group_id = fingerprint(record.source_samples)
        grouped.setdefault(record.source_group_id, []).append(record)

    groups = {}
    for group_id, members in grouped.items():
        eligible = [r for r in members if _is_eligible(r)]
        studies = sorted({r.study_accession for r in eligible if r.study_accession})
        # DECLARED runs, not the source sample's runs. The latter are identical for every
        # member of the group, so comparing them would make n_distinct_run_sets always 1
        # and the two run-mismatch verdicts unreachable.
        run_sets = {tuple(sorted(set(r.declared_runs))) for r in eligible if r.declared_runs}
        with_runs = sum(1 for r in members if r.declared_runs)
        group = SourceGroup(
            source_group_id=group_id,
            source_samples=sorted(set(members[0].source_samples)),
            n_genomes=len(members),
            n_events=len(studies),
            independent=len(studies) >= 2,
            # Only a group that could hold an independent duplicate is worth comparing
            # genomes for. One study cannot, however similar its genomes are.
            needs_organism_check=len(studies) >= 2,
            n_distinct_run_sets=len(run_sets),
            run_mismatch=len(run_sets) > 1,
            run_coverage=(float(with_runs) / len(members)) if members else 0.0,
            n_ineligible=len(members) - len(eligible),
            studies=studies,
            centres=sorted({r.center_name for r in members if r.center_name}),
        )
        groups[group_id] = group
        for record in members:
            record.relationship = _verdict(group, record)
    return groups


def _verdict(group, record):
    """The verdict from provenance alone. Never `duplicate` -- that needs the join."""
    if record.source_partial:
        return REL_PARTIAL_SOURCE
    if group.n_events >= 2:
        # A declared-run mismatch does not change this: two runs of one library are
        # technical replicates, so it is recorded as `run_mismatch` and nothing more.
        return REL_PENDING_ORGANISM_CHECK
    if group.n_genomes == 1:
        return REL_SINGLETON
    return REL_SAME_STUDY_REDUNDANCY


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
#
# Logging to stderr, data to files. Rows are sorted, so output is diffable between runs.

MAGS_COLUMNS = [
    "source_group_id", "genome_accession", "version", "genome_type", "mag_sample",
    "in_input",
    "source_samples", "n_source_samples", "source_parent_type",
    # what THIS submitter said they used, and how we know
    "declared_runs", "run_evidence",
    "study_accession", "center_name", "last_updated", "first_public",
    "assembly_software", "binning_software",
    # Submitter-produced
    "completeness_score_submitter", "contamination_score_submitter",
    "assembly_type",
    # evidence only -- never an input to any decision
    "tax_id", "scientific_name",
    "n_events_in_group", "independent",
    # a partial source set keeps its group but cannot create a generation event
    "source_partial",
    # organism identity, blank until the gemsparcl join
    "gcu_id", "gcu_note",
    "relationship", "resolution_status",
]

SOURCE_GROUPS_COLUMNS = [
    "source_group_id", "source_samples", "n_genomes", "n_events", "independent",
    "needs_organism_check", "n_ineligible", "run_mismatch", "n_distinct_run_sets",
    "run_coverage",
    "studies", "centres",
]


def _bool(value):
    return "TRUE" if value else "FALSE"


def _join(values):
    return ",".join(values)


def write_mags_tsv(path, records, source_groups=None):
    """One row per genome, INCLUDING every input that could not be resolved.

    An unresolved input keeps its row with a `resolution_status` saying why, so nothing
    disappears between the input file and the output.
    """
    source_groups = source_groups or {}
    ordered = sorted(records,
                     key=lambda r: (r.source_group_id, r.accession, r.input_accession))
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MAGS_COLUMNS, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for record in ordered:
            group = source_groups.get(record.source_group_id)
            writer.writerow({
                "source_group_id": record.source_group_id,
                "genome_accession": record.accession,
                "version": record.version,
                "genome_type": record.genome_type,
                "mag_sample": record.mag_sample,
                "in_input": _bool(record.in_input),
                "source_samples": _join(sorted(set(record.source_samples))),
                "n_source_samples": str(len(set(record.source_samples))),
                "source_parent_type": _join(_dedupe(record.source_parent_types)),
                "declared_runs": _join(sorted(set(record.declared_runs))),
                "run_evidence": record.run_evidence,
                "study_accession": record.study_accession,
                "center_name": record.center_name,
                "last_updated": record.last_updated,
                "first_public": record.first_public,
                "assembly_software": record.assembly_software,
                "binning_software": record.binning_software,
                "completeness_score_submitter": record.completeness_score,
                "contamination_score_submitter": record.contamination_score,
                "assembly_type": record.assembly_type,
                "tax_id": record.tax_id,
                "scientific_name": record.scientific_name,
                "n_events_in_group": str(group.n_events) if group else "",
                "independent": _bool(group.independent) if group else "",
                "source_partial": _bool(record.source_partial),
                "gcu_id": record.gcu_id,
                "gcu_note": record.gcu_note,
                "relationship": record.relationship,
                "resolution_status": record.resolution_status,
            })


def write_source_groups_tsv(path, source_groups):
    """One row per source group: genomes declaring the same source sample(s)."""
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SOURCE_GROUPS_COLUMNS, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for group_id in sorted(source_groups):
            group = source_groups[group_id]
            writer.writerow({
                "source_group_id": group.source_group_id,
                "source_samples": _join(group.source_samples),
                "n_genomes": str(group.n_genomes),
                "n_events": str(group.n_events),
                "independent": _bool(group.independent),
                "needs_organism_check": _bool(group.needs_organism_check),
                "n_ineligible": str(group.n_ineligible),
                "run_mismatch": _bool(group.run_mismatch),
                "n_distinct_run_sets": str(group.n_distinct_run_sets),
                "run_coverage": "%.3f" % group.run_coverage,
                "studies": _join(group.studies),
                "centres": _join(group.centres),
            })


def summarise(records, source_groups):
    """Run statistics. Resolution rate is the headline metric."""
    inputs = [r for r in records if r.in_input]
    resolved = [r for r in inputs if r.resolution_status == STATUS_RESOLVED]
    traceable = [r for r in inputs
                 if r.resolution_status not in (STATUS_NOT_A_MAG, STATUS_INVALID)]
    status_counts = {}
    for record in inputs:
        status_counts[record.resolution_status] = (
            status_counts.get(record.resolution_status, 0) + 1)
    return {
        "inputs": len(inputs),
        "resolved_inputs": len(resolved),
        # Counted over INPUTS only: discovered rows are found by definition, so including
        # them would flatter the rate and hide how much of the user's list was traceable.
        "resolution_rate": (float(len(resolved)) / len(inputs)) if inputs else 0.0,
        # Recall over the inputs we tried to trace. `not_a_mag` and `invalid_accession`
        # are correct answers about the input, not failures, so they are excluded. Not the
        # same as "the inputs that are MAGs": a blank assembly_type is accepted.
        "traceable_inputs": len(traceable),
        "traceable_resolution_rate": (
            (float(len(resolved)) / len(traceable)) if traceable else 0.0),
        "discovered": sum(1 for r in records if not r.in_input),
        "source_groups": len(source_groups),
        "independent_groups": sum(1 for g in source_groups.values() if g.independent),
        "groups_needing_organism_check": sum(
            1 for g in source_groups.values() if g.needs_organism_check),
        "ineligible_rows": sum(1 for r in records if not _is_eligible(r)),
        "status_counts": status_counts,
    }


def log_summary(summary):
    LOG.info("input accessions:        %s", summary["inputs"])
    LOG.info("resolved inputs:         %s", summary["resolved_inputs"])
    LOG.info("RESOLUTION RATE:         %.1f%% of all inputs", 100.0 * summary["resolution_rate"])
    if summary["traceable_inputs"] != summary["inputs"]:
        LOG.info("  of the %s we tried to trace: %.1f%%   <- the provenance recall figure",
                 summary["traceable_inputs"],
                 100.0 * summary["traceable_resolution_rate"])
    LOG.info("discovered (not input):  %s", summary["discovered"])
    LOG.info("source groups:           %s", summary["source_groups"])
    LOG.info("  independent:           %s", summary["independent_groups"])
    if summary["source_groups"]:
        # How often independent re-analysis actually happens -- the headline answer.
        LOG.info("  = %.1f%% of source groups hold genomes from more than one study",
                 100.0 * summary["independent_groups"] / float(summary["source_groups"]))
    LOG.info("needing an organism check: %s", summary["groups_needing_organism_check"])
    if summary["ineligible_rows"]:
        LOG.info("rows excluded from the independence call (partial source): %s",
                 summary["ineligible_rows"])
    for status, count in sorted(summary["status_counts"].items()):
        LOG.info("  status %-32s %s", status, count)
    LOG.info("NOTE: a source group is not a duplicate set. Organism identity is unknown "
             "until the gemsparcl join.")


# ---------------------------------------------------------------------------
# Tier 1: attach a gemsparcl cluster
# ---------------------------------------------------------------------------
#
# Where ORGANISM identity enters. The script never RUNS gemsparcl; it reads gemsparcl's
# files and lists whatever they could not answer. Two sources of a cluster, cheapest first:
#
#   wdp_clusters.csv   already computed over 3.28M ENA genomes -> a LOOKUP
#   gemsparcl query    for genomes absent from that run        -> incremental
#
# A cluster means SAME ORGANISM and nothing about provenance. NOTHING here may set a
# duplicate verdict.

ENA_FASTA_URL = "https://www.ebi.ac.uk/ena/browser/api/fasta/%s?download=true"

NOTE_NO_HIT = "no_hit"


@dataclass(frozen=True)
class GcuAssignment:
    """One genome's organism assignment, as gemsparcl reported it."""

    genome_id: str
    gcu_id: str
    note: str = ""


@dataclass(frozen=True)
class CandidateEntry:
    """One genome still needing an organism check, and where to fetch it."""

    genome_accession: str
    source_group_id: str
    fasta_url: str


def read_gemsparcl_clusters(path, *, wanted, namespace="wdp"):
    """Read gemsparcl's cluster assignments into {genome_accession: GcuAssignment}.

    Two layouts are accepted, because gemsparcl writes a different one in each mode:

        genome_id,cluster              from `gemsparcl cluster`
        GCA_036238915,1

        genome_id,GCU,note             from `gemsparcl query`
        GCA_036238915,1,assigned

    A file with neither layout is rejected with a ValueError.

    Args:
        path: a gemsparcl `*_clusters.csv` or `*_query_results.csv`.
        namespace: prefix for every cluster id, so "wdp" turns cluster 1 into "wdp:1".
            Needed because the ids are plain integers meaning "the Nth cluster of this
            run" -- cluster 1 in two different files is two different organisms, and
            joining them unprefixed would report a duplicate that does not exist.
        wanted: the genome accessions to keep. Rows for anything else are skipped as the
            file is read. Pass None to keep the whole file.

            This has no default on purpose. The file describes every genome WDP has
            clustered -- around 3.3 million rows -- while a run needs a few thousand.
            Reading it whole exhausted the memory on an HPC node and the job was killed
            without an error message, so it is better to have to say what you want.

    A row with an empty cluster cell gets an empty `gcu_id`. That is what `gemsparcl
    query` writes for a `no_hit`, and it is an answer rather than a gap: the genome
    matched no existing cluster, so it is a novel organism and nothing can duplicate it.
    """
    assignments = {}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        columns = {name.lower(): name for name in (reader.fieldnames or ())}
        id_col = columns.get("genome_id")
        gcu_col = columns.get("cluster") or columns.get("gcu")
        if not id_col or not gcu_col:
            raise ValueError(
                "%s: expected genome_id and cluster columns, found %s"
                % (path, reader.fieldnames))
        note_col = columns.get("note")
        for row in reader:
            genome = (row.get(id_col) or "").strip()
            if not genome or (wanted is not None and genome not in wanted):
                continue
            raw = (row.get(gcu_col) or "").strip()
            assignments[genome] = GcuAssignment(
                genome_id=genome,
                gcu_id=(namespace + ":" + raw) if raw else "",
                note=(row.get(note_col) or "").strip() if note_col else "",
            )
    return assignments


def attach_gcus(records, assignments):
    """Copy cluster ids onto records."""
    attached = 0
    for record in records:
        entry = assignments.get(record.accession)
        if entry is None:
            continue
        record.gcu_id = entry.gcu_id
        record.gcu_note = entry.note
        attached += 1
    return attached


def _fasta_url(record):
    """Where to fetch this genome's sequence."""
    if record.genome_type == KIND_ERZ and record.generated_ftp:
        ftp = record.generated_ftp.split(";")[0].strip()
        return ftp if ftp.startswith("http") else "https://" + ftp
    return ENA_FASTA_URL % record.accession


def build_organism_check_candidates(records, groups):
    """Genomes whose organism is still unknown, with a URL to fetch each one.

    Skips single-study groups, ineligible rows, and rows that already have a `gcu_id`.
    """
    candidates = []
    for record in sorted(records, key=lambda r: (r.source_group_id, r.accession)):
        group = groups.get(record.source_group_id)
        if group is None or not group.needs_organism_check or not _is_eligible(record):
            continue
        if record.gcu_id:
            continue
        candidates.append(CandidateEntry(record.accession, record.source_group_id,
                                         _fasta_url(record)))
    return candidates


CANDIDATES_FILE = "organism_check_candidates.tsv"
CANDIDATE_COLUMNS = ["source_group_id", "genome_accession", "fasta_url"]


def write_organism_check_candidates_tsv(path, candidates):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(CANDIDATE_COLUMNS)
        for entry in candidates:
            writer.writerow([entry.source_group_id, entry.genome_accession,
                             entry.fasta_url])


# ---------------------------------------------------------------------------
# Tier 2: join organism identity to provenance -- the final verdict
# ---------------------------------------------------------------------------
#
#   same source group  (same declared source sample, not necessarily the same runs)
#     AND same cluster (same organism)
#     AND >= 2 studies (different submitters)
#     = DUPLICATE

REL_DUPLICATE_INDEPENDENT = "duplicate_independent"
REL_DUPLICATE_SAME_STUDY = "duplicate_same_study"
REL_CO_DERIVED_DIFFERENT_ORGANISM = "co_derived_different_organism"
REL_ORGANISM_UNKNOWN = "organism_unknown"

DUPLICATES_COLUMNS = [
    "duplicate_id", "source_group_id", "gcu_id", "n_genomes", "n_studies",
    "independent", "min_internal_ani", "genomes", "studies", "centres", "gcu_notes",
]


@dataclass
class DuplicateGroup:
    """Same source sample AND same organism. Independent when >= 2 studies."""

    duplicate_id: str
    source_group_id: str
    gcu_id: str
    genomes: List[str]
    studies: List[str]
    centres: List[str]
    n_studies: int
    independent: bool
    min_internal_ani: Optional[float] = None
    gcu_notes: List[str] = field(default_factory=list)


def read_gemsparcl_dists(path, *, wanted):
    """Read how similar gemsparcl measured each pair of genomes to be.

    This is what lets the output say how alike two genomes in a duplicate group actually
    are, instead of only that they ended up in the same cluster.

    gemsparcl writes one tab-separated row per comparison:

        GCA_937896655   GCA_977328445   0.9977

    The score is average nucleotide identity, where 1.0 means identical and gemsparcl
    treats 0.98 and above as the same organism. The returned dictionary is keyed by the
    pair of accessions without regard to order, so a lookup works whichever way round the
    two genomes are given.

    Not every pair is in the file, and that is normal: gemsparcl only records each genome
    against its closest neighbours, 50 of them by default. A pair that is missing was not
    among those, which is not the same as being dissimilar -- so a missing pair reads as
    "not measured" and min_internal_ani returns None rather than 0.0. To get every pair
    inside one source group, re-run gemsparcl on that group alone.

    Args:
        path: a gemsparcl `.dists` file.
        wanted: the genome accessions to keep. Pass None to keep the whole file.
    """
    dists = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for raw in handle:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                value = float(parts[2])
            except ValueError:
                continue  # header, or a malformed line
            left, right = parts[0].strip(), parts[1].strip()
            if wanted is not None and (left not in wanted or right not in wanted):
                continue
            dists[frozenset((left, right))] = value
    return dists


def min_internal_ani(genomes, dists):
    """Lowest ANI between any two members, or None if no pair was measured.

    The minimum, so a loosely-connected group is visible: clustering links genomes in a
    chain, so A-B and B-C can both clear the threshold while A-C does not.
    """
    values = []
    for i in range(len(genomes)):
        for j in range(i + 1, len(genomes)):
            value = dists.get(frozenset((genomes[i], genomes[j])))
            if value is not None:
                values.append(value)
    return min(values) if values else None


def assign_duplicates(records, groups, dists=None):
    """Join organism identity onto provenance and set the final verdict."""
    dists = dists or {}
    buckets = {}
    for record in records:
        group = groups.get(record.source_group_id)
        if group is None or not group.needs_organism_check or not _is_eligible(record):
            continue
        if not record.gcu_id:
            # `no_hit` is an ANSWER: the organism matched nothing, so it is novel and
            # nothing can duplicate it. A blank with no note means we could not look.
            # Those are different claims and must not share a label.
            record.relationship = (REL_CO_DERIVED_DIFFERENT_ORGANISM
                                   if record.gcu_note == NOTE_NO_HIT
                                   else REL_ORGANISM_UNKNOWN)
            continue
        buckets.setdefault((record.source_group_id, record.gcu_id), []).append(record)

    duplicates = []
    for (group_id, gcu_id), members in sorted(buckets.items()):
        if len(members) < 2:
            # The only genome of this organism in the group. The rest came from the same
            # sample but are different organisms, so none of them duplicates it.
            for record in members:
                record.relationship = REL_CO_DERIVED_DIFFERENT_ORGANISM
            continue
        studies = sorted({r.study_accession for r in members if r.study_accession})
        independent = len(studies) >= 2
        verdict = REL_DUPLICATE_INDEPENDENT if independent else REL_DUPLICATE_SAME_STUDY
        for record in members:
            record.relationship = verdict
        accessions = sorted(r.accession for r in members)
        duplicates.append(DuplicateGroup(
            duplicate_id="dup_" + hashlib.sha1(
                "|".join(accessions).encode("utf-8")).hexdigest()[:12],
            source_group_id=group_id,
            gcu_id=gcu_id,
            genomes=accessions,
            studies=studies,
            centres=sorted({r.center_name for r in members if r.center_name}),
            n_studies=len(studies),
            independent=independent,
            min_internal_ani=min_internal_ani(accessions, dists),
            gcu_notes=sorted({r.gcu_note for r in members if r.gcu_note}),
        ))
    return duplicates


def write_duplicates_tsv(path, duplicates):
    """The answer. Independent duplicates first, since those are what was asked for."""
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(DUPLICATES_COLUMNS)
        for dup in sorted(duplicates, key=lambda d: (not d.independent, d.duplicate_id)):
            writer.writerow([
                dup.duplicate_id, dup.source_group_id, dup.gcu_id,
                str(len(dup.genomes)), str(dup.n_studies), _bool(dup.independent),
                "" if dup.min_internal_ani is None else "%.4f" % dup.min_internal_ani,
                _join(dup.genomes), _join(dup.studies), _join(dup.centres),
                ";".join(dup.gcu_notes),
            ])


# ---------------------------------------------------------------------------
# Reading our own output back (phase 2 picks up where phase 1 stopped)
# ---------------------------------------------------------------------------

def _unjoin(value):
    return [part for part in (value or "").split(",") if part]


def _unbool(value):
    return (value or "").strip().upper() == "TRUE"


def read_mags_tsv(path):
    """Rebuild GenomeRecords from a phase-1 mags.tsv."""
    records = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            record = GenomeRecord(
                accession=row.get("genome_accession", ""),
                input_accession=row.get("genome_accession", ""),
                genome_type=row.get("genome_type", ""),
                resolution_status=row.get("resolution_status", ""),
                version=row.get("version", ""),
                mag_sample=row.get("mag_sample", ""),
                study_accession=row.get("study_accession", ""),
                center_name=row.get("center_name", ""),
                last_updated=row.get("last_updated", ""),
                first_public=row.get("first_public", ""),
                assembly_software=row.get("assembly_software", ""),
                binning_software=row.get("binning_software", ""),
                completeness_score=row.get("completeness_score_submitter", ""),
                contamination_score=row.get("contamination_score_submitter", ""),
                assembly_type=row.get("assembly_type", ""),
                tax_id=row.get("tax_id", ""),
                scientific_name=row.get("scientific_name", ""),
            )
            record.in_input = _unbool(row.get("in_input"))
            record.source_samples = _unjoin(row.get("source_samples"))
            record.source_parent_types = _unjoin(row.get("source_parent_type"))
            record.declared_runs = _unjoin(row.get("declared_runs"))
            record.run_evidence = row.get("run_evidence", "")
            record.source_partial = _unbool(row.get("source_partial"))
            record.source_group_id = row.get("source_group_id", "")
            record.relationship = row.get("relationship", "")
            record.gcu_id = row.get("gcu_id", "")
            record.gcu_note = row.get("gcu_note", "")
            records.append(record)
    return records


def read_source_groups_tsv(path):
    """Rebuild SourceGroups from a phase-1 source_groups.tsv."""
    groups = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            group_id = row.get("source_group_id", "")
            if not group_id:
                continue
            groups[group_id] = SourceGroup(
                source_group_id=group_id,
                source_samples=_unjoin(row.get("source_samples")),
                n_genomes=int(row.get("n_genomes") or 0),
                n_events=int(row.get("n_events") or 0),
                independent=_unbool(row.get("independent")),
                needs_organism_check=_unbool(row.get("needs_organism_check")),
                n_distinct_run_sets=int(row.get("n_distinct_run_sets") or 0),
                run_coverage=float(row.get("run_coverage") or 0.0),
                run_mismatch=_unbool(row.get("run_mismatch")),
                n_ineligible=int(row.get("n_ineligible") or 0),
                studies=_unjoin(row.get("studies")),
                centres=_unjoin(row.get("centres")),
            )
    return groups


def read_gemsparcl_results(directory, *, wanted):
    """Read every gemsparcl output in a directory, namespaced by file.

    Each file's cluster ids are namespaced with its own stem, because the integers mean
    "the Nth cluster of THIS run": a per-group re-cluster and the global run both contain
    a cluster `1`, for different organisms.
    """
    assignments = {}
    dists = {}
    found = False
    for name in sorted(os.listdir(directory)):
        path = os.path.join(directory, name)
        if name.endswith("_clusters.csv") or name.endswith("_query_results.csv"):
            stem = name[:-len("_clusters.csv")] if name.endswith("_clusters.csv") \
                else name[:-len("_query_results.csv")]
            assignments.update(read_gemsparcl_clusters(path, namespace=stem or name,
                                                       wanted=wanted))
            found = True
        elif name.endswith(".dists"):
            dists.update(read_gemsparcl_dists(path, wanted=wanted))
            found = True
    if not found:
        raise ValueError(
            "%s contains no gemsparcl output (*_clusters.csv, *_query_results.csv, "
            "*.dists). Refusing to report 'no duplicates' for an empty directory."
            % directory)
    return assignments, dists


# ---------------------------------------------------------------------------
# CLI and orchestration
# ---------------------------------------------------------------------------

def read_accessions(path):
    """One accession per line; blanks, '#' comments and extra columns ignored.

    Takes an NCBI assembly_summary file unchanged. Order is preserved, duplicates dropped.
    """
    out, seen = [], set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.split("#", 1)[0]
            token = line.split("\t")[0].split(",")[0].strip()
            if token and token not in seen:
                seen.add(token)
                out.append(token)
    return out


def _resolve_discovered(client, discovered_samples):
    """Resolve discovered samples to genomes; those with none are `sample_only`.

    Asks `assembly` and `analysis`, both. For an input accession we know which of the two
    holds it, but a discovered sample's genome could be at either level, and asking only
    `assembly` would miss every cohort deposited as an ERZ analysis. One benchmark cohort
    is 46 genomes, all of them ERZ.

    Uses `query=sample_accession=...` rather than `includeAccessions`, which rejects a
    batch containing both accession forms -- and a batch of discovered samples routinely
    contains both.
    """
    samples = list(discovered_samples)
    by_sample = {}
    for start in range(0, len(samples), DEFAULT_OR_BATCH):
        chunk = samples[start:start + DEFAULT_OR_BATCH]
        query = " OR ".join('sample_accession="%s"' % s for s in chunk)
        for result, fields, acc_key, kind in (
                ("assembly", ASSEMBLY_FIELDS, "accession", KIND_GCA),
                ("analysis", ANALYSIS_FIELDS, "analysis_accession", KIND_ERZ)):
            for row in client.search(result, fields, query=query):
                record = _row_to_record(row, acc_key, row.get(acc_key, ""), kind)
                record.in_input = False
                by_sample.setdefault(record.mag_sample, []).append(record)

    out = []
    for sample in samples:
        candidates = by_sample.get(sample)
        if candidates:
            out.append(_latest(candidates, lambda r: r.version))
        else:
            out.append(GenomeRecord(accession="", input_accession="", genome_type="",
                                    resolution_status=STATUS_SAMPLE_ONLY,
                                    mag_sample=sample))
            out[-1].in_input = False
    return out


def _apply_resolution(record, resolution):
    """Copy a SourceResolution onto a record.

    `source_partial` is LOAD-BEARING: without it the partial-source guard is dead code and
    a co-assembly whose source set is a known subset of the truth looks complete, so it
    groups alone and takes its generation event with it.
    """
    if resolution is None:
        return
    if resolution.status == STATUS_RESOLVED:
        record.source_samples = list(resolution.source_samples)
        record.source_parent_types = list(resolution.parent_types)
        record.source_sample_runs = list(resolution.runs)
        record.source_partial = resolution.partial
    elif record.resolution_status == STATUS_RESOLVED:
        record.resolution_status = resolution.status


def run_pipeline(client, accessions, expand_enabled=True):
    """Phase 1: resolve, walk up, expand, resolve discovered, group.

    Returns (records, source_groups). A per-accession failure becomes a `resolution_status`
    value and the run completes with counts, never a mid-loop exit.
    """
    records = resolve_genomes(client, accessions)
    inputs = [records[a.raw] for a in accessions]
    mags = [r for r in inputs if r.resolution_status == STATUS_RESOLVED and r.mag_sample]

    sources = resolve_sources(client, [r.mag_sample for r in mags]) if mags else {}
    for record in mags:
        _apply_resolution(record, sources.get(record.mag_sample))

    all_records = list(inputs)
    if expand_enabled and sources:
        keys = expansion_keys(sources.values())
        known = {r.mag_sample for r in mags}
        discovered = [s for s in expand(client, keys) if s not in known]
        if discovered:
            found = _resolve_discovered(client, discovered)
            resolved = [r for r in found if r.resolution_status == STATUS_RESOLVED]
            extra_sources = (resolve_sources(client, [r.mag_sample for r in resolved])
                             if resolved else {})
            for record in resolved:
                _apply_resolution(record, extra_sources.get(record.mag_sample))
            all_records += found

    return all_records, build_source_groups(all_records)


def _phase1(args, client):
    accessions = [classify_accession(a) for a in read_accessions(args.accessions)]
    LOG.info("read %s accession(s) from %s", len(accessions), args.accessions)
    if not accessions:
        LOG.error("no accessions found in %s", args.accessions)
        return 1

    records, groups = run_pipeline(client, accessions,
                                   expand_enabled=not args.no_expand)

    known = {}
    if args.gemsparcl_clusters:
        # Discovered genomes need their clusters too, not just the input ones -- most
        # duplicate pairs are between two discovered genomes.
        known = read_gemsparcl_clusters(
            args.gemsparcl_clusters,
            wanted={r.accession for r in records if r.accession})
        LOG.info("tier 1: %s genome(s) already have a cluster from %s",
                 attach_gcus(records, known), args.gemsparcl_clusters)

    candidates = build_organism_check_candidates(records, groups)
    write_mags_tsv(os.path.join(args.outdir, "mags.tsv"), records, groups)
    write_source_groups_tsv(os.path.join(args.outdir, "source_groups.tsv"), groups)
    write_organism_check_candidates_tsv(
        os.path.join(args.outdir, CANDIDATES_FILE), candidates)

    summary = summarise(records, groups)
    log_summary(summary)
    LOG.info("wrote mags.tsv, source_groups.tsv, %s to %s", CANDIDATES_FILE, args.outdir)
    if candidates:
        LOG.info("%s genome(s) have no cluster yet; they are listed in %s and will "
                 "report 'organism_unknown'", len(candidates), CANDIDATES_FILE)
    return 0


def _phase2(args):
    mags_path = os.path.join(args.outdir, "mags.tsv")
    groups_path = os.path.join(args.outdir, "source_groups.tsv")
    for path in (mags_path, groups_path):
        if not os.path.exists(path):
            # "No duplicates found" and "nothing to look at" must not share an exit path.
            LOG.error("%s not found. Run phase 1 first (--accessions).", path)
            return 1

    records = read_mags_tsv(mags_path)
    groups = read_source_groups_tsv(groups_path)
    wanted = {r.accession for r in records if r.accession}
    assignments, dists = read_gemsparcl_results(args.gemsparcl_results, wanted=wanted)
    LOG.info("read %s cluster assignment(s) and %s pairwise distance(s)",
             len(assignments), len(dists))
    LOG.info("tier 1: attached a cluster to %s genome(s)",
             attach_gcus(records, assignments))

    duplicates = assign_duplicates(records, groups, dists=dists)
    write_duplicates_tsv(os.path.join(args.outdir, "duplicates.tsv"), duplicates)
    write_mags_tsv(mags_path, records, groups)
    # Written twice on purpose. Phase 1 knew no clusters yet, so it listed every
    # candidate; after the join, only the genomes still without one are left.
    write_organism_check_candidates_tsv(
        os.path.join(args.outdir, CANDIDATES_FILE),
        build_organism_check_candidates(records, groups))

    independent = [d for d in duplicates if d.independent]
    LOG.info("DUPLICATE GROUPS:        %s", len(duplicates))
    LOG.info("  independent:           %s", len(independent))
    LOG.info("  same-study redundancy: %s", len(duplicates) - len(independent))
    LOG.info("  genomes involved:      %s",
             sum(len(d.genomes) for d in independent))
    unknown = sum(1 for r in records if r.relationship == REL_ORGANISM_UNKNOWN)
    if unknown:
        LOG.warning("%s genome(s) had NO cluster assignment -- 'organism_unknown' is not "
                    "the same as 'not a duplicate'", unknown)
    chained = [d for d in independent
               if d.min_internal_ani is not None and d.min_internal_ani < 0.99]
    if chained:
        LOG.warning("%s duplicate group(s) are chained (an internal pair below 0.99) -- "
                    "see min_internal_ani", len(chained))
    LOG.info("wrote duplicates.tsv to %s", args.outdir)
    return 0


def build_arg_parser():
    parser = argparse.ArgumentParser(
        prog="find_mag_source_duplicates.py",
        description="Find MAGs that are the same organism, built from the same source "
                    "by different groups. Runs in two phases with a gemsparcl step "
                    "between them.")
    parser.add_argument("--accessions",
                        help="file of INSDC genome accessions, one per line")
    parser.add_argument("--gemsparcl-results",
                        help="directory of gemsparcl output (*_clusters.csv, *.dists) to "
                             "check candidates against")
    parser.add_argument("--outdir", default="results",
                        help="where the output files live (default: results)")
    parser.add_argument("--gemsparcl-clusters",
                        help="phase 1, optional: an existing wdp_clusters.csv to look up")
    parser.add_argument("--cache-dir",
                        help="on-disk response cache; the HPC resumability mechanism")
    parser.add_argument("--no-cache", action="store_true", help="disable the cache")
    parser.add_argument("--no-expand", action="store_true",
                        help="do not look beyond the input list (much lower recall)")
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--log-level", default="INFO")
    return parser


def main(argv=None):
    args = build_arg_parser().parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO),
                        format="%(levelname)s %(message)s", stream=sys.stderr)
    if not args.accessions and not args.gemsparcl_results:
        LOG.error("give --accessions (phase 1) or --gemsparcl-results (phase 2)")
        return 2
    try:
        os.makedirs(args.outdir, exist_ok=True)
        if args.accessions:
            cache = None if args.no_cache else args.cache_dir
            client = EnaClient(cache_dir=cache, timeout=args.timeout,
                               batch_size=args.batch_size)
            code = _phase1(args, client)
            if code or not args.gemsparcl_results:
                return code
        return _phase2(args)
    except (EnaError, ValueError, OSError) as exc:
        LOG.error("%s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
