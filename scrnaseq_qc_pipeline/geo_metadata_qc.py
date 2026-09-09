"""Cross-check curated Synapse annotations against the sample's own GEO record.

Schema validation at curation time only confirms a value is *legal* (drawn
from the controlled vocabulary) — it can't confirm a curator picked the
*correct* one for this specific GEO-derived sample (e.g. nothing stops
"Illumina HiSeq" from being entered when GEO's own SOFT record says
NovaSeq). GEO is an independent, authoritative provenance source for a
handful of fields, so this module flags divergence between what's on
Synapse and what GEO itself says for the same GSM.

Only fields GEO actually records independently are checked here — platform,
reference genome, and library-prep chemistry. Fields like DataLevel or
Component are ADA/Synapse bookkeeping with no GEO equivalent and are out of
scope for this comparison.
"""

from __future__ import annotations

from typing import Dict, List, Optional

from qc_common import Finding, first_value, make_finding

try:
    from geo_synapse.geo import combine_metadata_frames, extract_geo_metadata
    from geo_synapse.ena import map_platform
    HAS_GEO_SYNAPSE = True
except ImportError:
    HAS_GEO_SYNAPSE = False

GEO_SYNAPSE_INSTALL_HINT = (
    "geo-synapse is required to fetch GEO sample metadata. "
    "pip install git+https://github.com/sagebio-ada/geo_dataset_creation.git"
)


def fetch_geo_metadata_by_gsm(geo_accession: str) -> Dict[str, dict]:
    """Return {GSM id: row dict} for every sample in a GEO series."""
    frames = extract_geo_metadata(geo_accession)
    df = combine_metadata_frames(frames)
    if df.empty or "gsm_id" not in df.columns:
        return {}
    return {str(row["gsm_id"]): row.to_dict() for _, row in df.iterrows()}


def _norm(s: object) -> str:
    return str(s or "").strip().lower()


def _norm_compact(s: object) -> str:
    return _norm(s).replace(" ", "").replace("-", "").replace("_", "")


def _get_kit_text(geo_row: dict) -> str:
    """The reagent-kit column's apostrophe varies (curly vs straight) by feed."""
    for key, val in geo_row.items():
        if "reagent_kit" in key.lower():
            return str(val or "")
    return ""


def check_geo_annotations(
    entity_id: str, entity_name: str, annotations: Dict[str, list], geo_row: dict
) -> List[Finding]:
    findings: List[Finding] = []

    def add(check: str, status: str, detail: str) -> None:
        findings.append(make_finding(entity_id, entity_name, check, status, detail))

    # No pass/fail verdict on any of these -- whether a string/substring match
    # counts as "correct" is a judgment call, and these are our own comparison,
    # not a tool's. Report both sides; let a human decide if they agree.
    declared_platform = first_value(annotations.get("platform"))
    geo_instrument = str(geo_row.get("sample_instrument_model") or "")
    if declared_platform and geo_instrument:
        mapped = map_platform(geo_instrument) or geo_instrument
        add("platform_vs_geo", "INFO",
            f"Synapse platform={declared_platform!r} vs GEO sample_instrument_model={geo_instrument!r} "
            f"(normalized: {_norm_compact(declared_platform)!r} vs {_norm_compact(mapped)!r})")

    declared_ref = first_value(annotations.get("referenceSet"))
    geo_assembly = str(geo_row.get("assembly") or "")
    if declared_ref and geo_assembly:
        add("reference_vs_geo", "INFO",
            f"Synapse referenceSet={declared_ref!r} vs GEO assembly={geo_assembly!r}")

    kit_text = _get_kit_text(geo_row)
    declared_version = first_value(annotations.get("libraryVersion"))
    if declared_version and kit_text:
        add("library_version_vs_geo", "INFO",
            f"Synapse libraryVersion={declared_version!r} vs GEO kit description={kit_text!r}")

    declared_method = first_value(annotations.get("libraryPreparationMethod"))
    if declared_method and kit_text:
        add("library_prep_method_vs_geo", "INFO",
            f"Synapse libraryPreparationMethod={declared_method!r} vs GEO kit description={kit_text!r}")

    return findings
