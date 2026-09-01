"""
compound_names.py

Shared helpers for deciding what a compound synonym actually is.

PubChem returns synonym lists that mix a handful of real names in with database
accessions, registry numbers, InChIKeys and full IUPAC strings. Two stages need
to tell those apart: picking a display name, and deciding which synonyms are
worth sending to ClinicalTrials.gov.
"""

import re

# Accession prefixes that only ever appear inside a database, never in a paper
# or a trial registration.
_ACCESSION_PREFIXES = (
    "SCHEMBL", "CHEMBL", "BDBM", "DTXSID", "DTXCID", "NCGC", "MFCD",
    "AC1", "CID", "SID", "REFCHEM", "ZINC", "AKOS", "CS-", "HY-",
    "DB", "NSC", "EINECS", "UNII", "US", "WO", "EP", "Q", "EN300", "BRD-",
    "CHEBI", "KBIO", "CBIOL", "GLXC", "EX-A", "BSPBIO", "SPECTRUM", "SMR",
    "MLS", "HMS", "STK", "IDI1", "NCIMECH",
)

# Patent numbers and vendor catalogue entries masquerading as names.
_JUNK_MARKERS = ("code no.", "patent", "cas-", "einecs", "component of")

_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")
_UNII_RE = re.compile(r"^[A-Z0-9]{10}$")
_ATC_RE = re.compile(r"^[A-Z]\d{2}[A-Z]{2}\d{2}$")
_NUMERIC_RE = re.compile(r"^[\d\s.,-]+$")

# Trailing source tags PubChem appends, e.g. "Trifarotene [INN]".
_TAG_RE = re.compile(r"\s*\[[^\]]*\]\s*$")

# Official-name tags. The bare language-neutral ones are the English INN.
_OFFICIAL_TAGS = ("[INN]", "[USAN:INN]", "[USAN]", "[BAN]", "[JAN]", "[INN:BAN]")

# Markers of a systematic/IUPAC name rather than a common one.
_IUPAC_MARKERS = ("[", "]", "{", "}")
_IUPAC_WORDS = (
    "carboxylic acid", "-yl)", "-yl]", "benzoic acid,", "amino-", "-oxy",
    "tetrahydro", "dimethylethyl", "nonatetraenoic", "cyclohex",
)


def strip_source_tag(name: str) -> str:
    """Remove a trailing '[INN]'-style tag, leaving the name itself."""
    prev = None
    out = (name or "").strip()
    while out != prev:
        prev = out
        out = _TAG_RE.sub("", out).strip()
    return out.strip(" ?,;:.-")


def official_name(candidates):
    """
    The INN/USAN name if PubChem tagged one, else None.

    Only untranslated tags count: '[INN-French]' gives the French spelling, which
    is not what we want to display or search.
    """
    for raw in candidates or []:
        s = (raw or "").strip()
        upper = s.upper()
        for tag in _OFFICIAL_TAGS:
            if upper.endswith(tag):
                name = strip_source_tag(s)
                if name and not is_accession(name):
                    return name
    return None


def is_accession(name: str) -> bool:
    """True for database identifiers, registry numbers and InChIKeys."""
    s = (name or "").strip()
    if not s:
        return True
    upper = s.upper()

    if _INCHIKEY_RE.match(upper) or _CAS_RE.match(s) or _ATC_RE.match(upper):
        return True
    if _NUMERIC_RE.match(s):
        return True
    # A bare 10-character alphanumeric block with digits is a UNII, but a real
    # word of that length is not.
    if _UNII_RE.match(upper) and any(c.isdigit() for c in upper) and any(c.isalpha() for c in upper):
        if not re.search(r"[aeiou]{2}|[A-Z][a-z]{3}", s):
            return True

    lowered = s.lower()
    if any(j in lowered for j in _JUNK_MARKERS):
        return True

    for prefix in _ACCESSION_PREFIXES:
        if upper.startswith(prefix):
            rest = upper[len(prefix):].lstrip(":-")
            # SCHEMBL381691 is an accession; "CS-Retinoid" would not be.
            if rest.isdigit() or (rest and rest.isalnum() and any(c.isdigit() for c in rest)):
                return True
    return False


def is_systematic(name: str) -> bool:
    """True for IUPAC / systematic chemical names."""
    s = (name or "").strip()
    if not s:
        return True
    if any(m in s for m in _IUPAC_MARKERS):
        return True
    lowered = s.lower()
    if any(w in lowered for w in _IUPAC_WORDS):
        return True
    # Heavily punctuated, digit-dense strings are systematic names.
    digits = sum(c.isdigit() for c in s)
    if len(s) > 45 and digits >= 3:
        return True
    if s.count(",") >= 3 and digits >= 3:
        return True
    return False


def is_searchable_term(name: str) -> bool:
    """
    Worth sending to ClinicalTrials.gov.

    Development codes such as CD5789 or BMS-189961 are kept: trials really are
    registered under them. Database accessions and IUPAC strings are not.
    """
    s = strip_source_tag(name)
    if len(s) < 3 or len(s) > 60:
        return False
    return not is_accession(s) and not is_systematic(s)


def _alpha_ratio(name):
    letters = sum(c.isalpha() for c in name)
    return letters / len(name) if name else 0.0


def is_bare_code(name):
    """A catalogue/development code such as 'E98923' rather than a word."""
    return _alpha_ratio(name) < 0.5 and any(c.isdigit() for c in name)


def expand_candidates(candidates):
    """
    Flatten synonym entries. PubChem sometimes packs several names into one
    string, e.g. 'Ro 13-7410;Arotinoid acid;AGN191183'.
    """
    out = []
    for raw in candidates or []:
        if not isinstance(raw, str):
            continue
        for part in raw.split(";"):
            part = part.strip()
            if part:
                out.append(part)
    return out


def choose_display_name(candidates):
    """
    Pick the most human-readable name from a list of synonyms, or None when the
    compound genuinely has no common name and is only known by accession.

    Returning None is deliberate: labelling a compound 'SCHEMBL382466' reads as a
    name when it is not one.
    """
    seen = set()
    usable = []
    for raw in expand_candidates(candidates):
        name = strip_source_tag(raw)
        if not name:
            continue
        key = name.lower()
        if key in seen:
            continue
        seen.add(key)
        if is_accession(name) or is_systematic(name):
            continue
        usable.append(name)

    if not usable:
        return None

    def rank(name):
        alpha = sum(c.isalpha() for c in name)
        digits = sum(c.isdigit() for c in name)
        # Real words first. Among the rest, favour the most alphabetic candidate:
        # 'BMS 189961' is a name, 'E98923' is a catalogue number, and ranking on
        # digit count alone would pick the catalogue number.
        return (
            0 if alpha >= max(3, len(name) - 2) else 1,
            -round(_alpha_ratio(name), 3),
            digits,
            len(name),
            name.lower(),
        )

    return sorted(usable, key=rank)[0]
