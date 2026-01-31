# rules/cross_checks.py
from typing import Any, Dict, List, Optional, Tuple
import os
import re

try:
    import pysam
    _HAVE_PYSAM = True
except Exception:
    pysam = None
    _HAVE_PYSAM = False


class Issue(object):
    def __init__(self, level, code, message, file, where=None, fix=None, verify=None):
        self.level = level  # "error" or "warning"
        self.code = code
        self.message = message
        self.file = file  # "cross"
        self.where = where
        self.fix = fix
        self.verify = verify

    def to_dict(self):
        return {
            "level": self.level,
            "code": self.code,
            "message": self.message,
            "file": self.file,
            "where": self.where,
            "fix": self.fix,
            "verify": self.verify,
        }


def _normalize_contig(name: str) -> str:
    """
    Used for *detecting* naming-style mismatch, not for rewriting.
    Examples:
      chr1 -> 1
      ChrM -> M
      chrMT -> MT
    """
    if name is None:
        return ""
    n = name.strip()
    n2 = re.sub(r"^(chr|CHR|Chr)", "", n)
    return n2


def _guess_naming_style(contigs: List[str]) -> str:
    """
    Returns: "chr" or "no_chr" or "mixed" or "unknown"
    """
    if not contigs:
        return "unknown"
    chr_like = 0
    no_chr_like = 0
    for c in contigs:
        if c.lower().startswith("chr"):
            chr_like += 1
        else:
            no_chr_like += 1
    if chr_like and no_chr_like:
        return "mixed"
    if chr_like:
        return "chr"
    return "no_chr"


def _finalize(issues: List[Issue], stats: Dict[str, Any]) -> Dict[str, Any]:
    errors = [i.to_dict() for i in issues if i.level == "error"]
    warnings = [i.to_dict() for i in issues if i.level == "warning"]
    return {"errors": errors, "warnings": warnings, "stats": stats}


def cross_check_fasta_vcf(
    fasta_contig_lengths: Dict[str, int],
    vcf_path: str,
    *,
    ref_check: bool = False,
    fasta_path: Optional[str] = None,
    max_records_scan: int = 200000,
    max_ref_checks: int = 2000,
) -> Dict[str, Any]:
    """
    Cross-check rules:
      1) VCF contig exists in FASTA
      2) VCF POS <= FASTA contig length
      3) Naming-style mismatch detection (chr1 vs 1)
      4) Optional: REF base matches FASTA base at POS (sampling)

    Params:
      fasta_contig_lengths: dict {contig_name: length}
      vcf_path: uploaded VCF path
      ref_check: if True, do sampled REF checks (requires fasta_path with index support)
      fasta_path: path to FASTA (required if ref_check=True)
      max_records_scan: performance limit on VCF record scanning
      max_ref_checks: limit on how many REF validations to do

    Returns:
      {"errors": [...], "warnings": [...], "stats": {...}}
    """
    issues: List[Issue] = []
    stats: Dict[str, Any] = {
        "vcf_records_scanned": 0,
        "truncated_scan": False,
        "missing_contig_count": 0,
        "out_of_range_pos_count": 0,
        "naming_mismatch_suspected": False,
        "fasta_style": None,
        "vcf_style": None,
        "ref_checks_enabled": bool(ref_check),
        "ref_checks_attempted": 0,
        "ref_mismatch_count": 0,
        "ref_check_skipped_reason": None,
        "parser_used": None,
        "missing_contigs": {},   # {contig: count}
        "out_of_range": {},      # {contig: count}
    }

    # Basic validation
    if not isinstance(fasta_contig_lengths, dict) or not fasta_contig_lengths:
        issues.append(Issue(
            level="error",
            code="FASTA_CONTIGS_MISSING",
            message="FASTA contig lengths are missing. Cross-check requires FASTA contig names and lengths.",
            file="cross",
            fix="Ensure your FASTA checks return a dict of contig lengths (e.g., {contig: length})."
        ))
        return _finalize(issues, stats)

    if not os.path.exists(vcf_path):
        issues.append(Issue(
            level="error",
            code="VCF_NOT_FOUND",
            message="VCF file path does not exist for cross-check.",
            file="cross",
            where=vcf_path
        ))
        return _finalize(issues, stats)

    fasta_contigs = list(fasta_contig_lengths.keys())
    stats["fasta_style"] = _guess_naming_style(fasta_contigs)

    # Optional REF check preparation
    fasta_fa = None
    if ref_check:
        if not fasta_path or not os.path.exists(fasta_path):
            stats["ref_checks_enabled"] = False
            stats["ref_check_skipped_reason"] = "FASTA path missing for REF check."
            issues.append(Issue(
                level="warning",
                code="REF_CHECK_SKIPPED_NO_FASTA",
                message="REF base check was requested but FASTA path was not provided; skipping REF validation.",
                file="cross",
                fix="Provide fasta_path when calling cross_check_fasta_vcf(ref_check=True, fasta_path=...)."
            ))
        elif not _HAVE_PYSAM:
            stats["ref_checks_enabled"] = False
            stats["ref_check_skipped_reason"] = "pysam not available for indexed FASTA access."
            issues.append(Issue(
                level="warning",
                code="REF_CHECK_SKIPPED_NO_PYSAM",
                message="REF base check skipped because pysam is not installed (needed for indexed FASTA access).",
                file="cross",
                fix="Install pysam and ensure FASTA is indexed (samtools faidx ref.fasta)."
            ))
        else:
            try:
                fasta_fa = pysam.FastaFile(fasta_path)
            except Exception as e:
                stats["ref_checks_enabled"] = False
                stats["ref_check_skipped_reason"] = "Could not open FASTA with pysam: {}".format(str(e))
                issues.append(Issue(
                    level="warning",
                    code="REF_CHECK_SKIPPED_FASTA_OPEN_FAIL",
                    message="REF base check skipped because FASTA could not be opened for random access.",
                    file="cross",
                    fix="Index the FASTA (samtools faidx) and ensure fasta_path points to the same file used to generate the VCF."
                ))

    # Run cross-check scan
    if _HAVE_PYSAM:
        stats["parser_used"] = "pysam"
        try:
            _scan_vcf_with_pysam(
                vcf_path,
                fasta_contig_lengths,
                issues,
                stats,
                ref_checks_enabled=bool(stats["ref_checks_enabled"]),
                fasta_fa=fasta_fa,
                max_records_scan=max_records_scan,
                max_ref_checks=max_ref_checks,
            )
        except Exception as e:
            issues.append(Issue(
                level="error",
                code="CROSSCHECK_VCF_PARSE_FAILED",
                message="Cross-check failed while reading VCF with pysam: {}".format(str(e)),
                file="cross",
                fix="Ensure the VCF is valid and readable. If it is bgzipped (.vcf.gz), ensure it is not corrupted."
            ))
            # fallback
            stats["parser_used"] = "manual_fallback"
            _scan_vcf_manually(
                vcf_path,
                fasta_contig_lengths,
                issues,
                stats,
                ref_checks_enabled=False,  # manual path doesn't do REF checks
                max_records_scan=max_records_scan,
            )
    else:
        stats["parser_used"] = "manual"
        _scan_vcf_manually(
            vcf_path,
            fasta_contig_lengths,
            issues,
            stats,
            ref_checks_enabled=False,
            max_records_scan=max_records_scan,
        )

    # Naming mismatch detection (heuristic)
    # If many missing contigs but their normalized forms match FASTA normalized forms -> likely chr mismatch.
    vcf_contigs_seen = list(stats.get("_vcf_contigs_seen", []))
    stats.pop("_vcf_contigs_seen", None)

    stats["vcf_style"] = _guess_naming_style(vcf_contigs_seen)

    if stats["missing_contig_count"] > 0 and fasta_contigs and vcf_contigs_seen:
        fasta_norm = set(_normalize_contig(c) for c in fasta_contigs)
        vcf_missing = list(stats["missing_contigs"].keys())
        missing_norm_match = 0
        for c in vcf_missing:
            if _normalize_contig(c) in fasta_norm:
                missing_norm_match += 1

        # If a significant fraction suggests normalization match, raise naming warning
        if missing_norm_match >= max(1, int(0.5 * len(vcf_missing))):
            stats["naming_mismatch_suspected"] = True
            issues.append(Issue(
                level="warning",
                code="CONTIG_NAMING_MISMATCH",
                message="VCF contig names appear to differ from FASTA naming style (e.g., chr1 vs 1).",
                file="cross",
                fix="Use consistent reference naming. Either rename VCF contigs to match FASTA or use the matching reference FASTA used for alignment/calling.",
                verify="Compare FASTA headers and VCF CHROM values; they should match exactly."
            ))

    # Summarize out-of-range / missing contigs as top-level issues
    if stats["missing_contig_count"] > 0:
        top = sorted(stats["missing_contigs"].items(), key=lambda x: x[1], reverse=True)[:5]
        issues.append(Issue(
            level="error",
            code="VCF_CONTIG_NOT_IN_FASTA",
            message="VCF contains contig(s) not present in FASTA. Top missing: {}".format(top),
            file="cross",
            fix="Use the correct reference FASTA for this VCF or rename contigs to match the FASTA headers."
        ))

    if stats["out_of_range_pos_count"] > 0:
        top = sorted(stats["out_of_range"].items(), key=lambda x: x[1], reverse=True)[:5]
        issues.append(Issue(
            level="error",
            code="VCF_POS_OUT_OF_RANGE",
            message="VCF contains position(s) beyond FASTA contig length. Top affected: {}".format(top),
            file="cross",
            fix="Check whether the VCF was called against a different reference build (e.g., GRCh37 vs GRCh38) or contig lengths differ."
        ))

    if stats.get("ref_checks_enabled") and stats["ref_mismatch_count"] > 0:
        issues.append(Issue(
            level="warning",
            code="REF_BASE_MISMATCH",
            message="Sampled REF base checks found {} mismatches out of {} attempts.".format(
                stats["ref_mismatch_count"], stats["ref_checks_attempted"]
            ),
            file="cross",
            fix="This often indicates a reference mismatch (different FASTA build) or left-alignment/normalization differences for indels.",
            verify="Confirm the VCF was generated against the same exact FASTA (same contig names and sequences)."
        ))

    return _finalize(issues, stats)


def _scan_vcf_with_pysam(
    vcf_path: str,
    fasta_contig_lengths: Dict[str, int],
    issues: List[Issue],
    stats: Dict[str, Any],
    *,
    ref_checks_enabled: bool,
    fasta_fa: Optional["pysam.FastaFile"],
    max_records_scan: int,
    max_ref_checks: int,
) -> None:
    v = pysam.VariantFile(vcf_path)
    vcf_contigs_seen = set()

    scanned = 0
    ref_attempts = 0
    ref_mismatch = 0

    for rec in v.fetch():
        scanned += 1
        if scanned > max_records_scan:
            stats["truncated_scan"] = True
            break

        chrom = str(rec.contig) if rec.contig is not None else ""
        pos = int(rec.pos) if rec.pos is not None else -1
        ref = str(rec.ref) if rec.ref is not None else ""

        if chrom:
            vcf_contigs_seen.add(chrom)

        # 1) contig exists
        if chrom not in fasta_contig_lengths:
            stats["missing_contig_count"] += 1
            stats["missing_contigs"][chrom] = stats["missing_contigs"].get(chrom, 0) + 1
            # no point checking range/ref if contig missing
            continue

        # 2) position within contig length
        clen = fasta_contig_lengths.get(chrom, None)
        if clen is not None and pos > clen:
            stats["out_of_range_pos_count"] += 1
            stats["out_of_range"][chrom] = stats["out_of_range"].get(chrom, 0) + 1
            continue

        # 3) optional REF base check (sampling, only for SNP-like REF length 1)
        if ref_checks_enabled and fasta_fa is not None and ref and len(ref) == 1 and pos >= 1:
            if ref_attempts >= max_ref_checks:
                continue
            ref_attempts += 1
            try:
                # pysam FastaFile uses 0-based, half-open intervals
                base = fasta_fa.fetch(chrom, pos - 1, pos).upper()
                if base and base != ref.upper():
                    ref_mismatch += 1
            except Exception:
                # Don't hard-fail; just skip
                pass

    stats["vcf_records_scanned"] = scanned if scanned <= max_records_scan else max_records_scan
    stats["ref_checks_attempted"] = ref_attempts
    stats["ref_mismatch_count"] = ref_mismatch
    stats["_vcf_contigs_seen"] = list(vcf_contigs_seen)


def _scan_vcf_manually(
    vcf_path: str,
    fasta_contig_lengths: Dict[str, int],
    issues: List[Issue],
    stats: Dict[str, Any],
    *,
    ref_checks_enabled: bool,
    max_records_scan: int,
) -> None:
    # Manual scan: only CHROM and POS from columns 1 and 2 after #CHROM.
    header_seen = False
    vcf_contigs_seen = set()

    scanned = 0
    with open(vcf_path, "r", encoding="utf-8", errors="replace") as f:
        for lineno, line in enumerate(f, start=1):
            line = line.rstrip("\n\r")
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_seen = True
                continue
            if not header_seen:
                continue

            scanned += 1
            if scanned > max_records_scan:
                stats["truncated_scan"] = True
                break

            parts = line.split("\t")
            if len(parts) < 2:
                continue

            chrom = parts[0].strip()
            pos_raw = parts[1].strip()

            if chrom:
                vcf_contigs_seen.add(chrom)

            # POS integer?
            try:
                pos = int(pos_raw)
            except Exception:
                # parsing issue belongs to VCF rule-checks; don't duplicate
                continue

            # contig exists
            if chrom not in fasta_contig_lengths:
                stats["missing_contig_count"] += 1
                stats["missing_contigs"][chrom] = stats["missing_contigs"].get(chrom, 0) + 1
                continue

            # range
            clen = fasta_contig_lengths.get(chrom, None)
            if clen is not None and pos > clen:
                stats["out_of_range_pos_count"] += 1
                stats["out_of_range"][chrom] = stats["out_of_range"].get(chrom, 0) + 1

    stats["vcf_records_scanned"] = scanned if scanned <= max_records_scan else max_records_scan
    stats["_vcf_contigs_seen"] = list(vcf_contigs_seen)
