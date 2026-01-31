# rules/vcf_checks.py
from typing import Any, Dict, List, Optional, Tuple
import os
import re

try:
    import pysam
    _HAVE_PYSAM = True
except Exception:
    pysam = None
    _HAVE_PYSAM = False


VALID_BASES = set("ACGTN")  # VCF REF/ALT commonly allow N. (Some pipelines allow other IUPAC codes; treat as warning.)


class Issue(object):
    def __init__(self, level, code, message, file, where=None, fix=None, verify=None):
        self.level = level  # "error" or "warning"
        self.code = code
        self.message = message
        self.file = file  # "vcf"
        self.where = where  # e.g. "line 123", "record chr1:100"
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


def _is_int(s: str) -> bool:
    try:
        int(s)
        return True
    except Exception:
        return False


def _is_float(s: str) -> bool:
    try:
        float(s)
        return True
    except Exception:
        return False


def _looks_like_iupac(seq: str) -> bool:
    # Basic check: contains letters beyond ACGTN. (Treat as warning, not always invalid.)
    return bool(re.search(r"[^ACGTN]", seq.upper()))


def _variant_type(ref: str, alt: str) -> str:
    # Very simple classification:
    # SNP: len(ref)==1 and len(alt)==1
    # INDEL: otherwise (includes MNP and complex; count as INDEL-like for your summary)
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    return "INDEL"


def vcf_rule_checks(vcf_path: str, max_records_scan: int = 200000) -> Dict[str, Any]:
    """
    Returns:
    {
        "errors": [Issue-as-dict...],
        "warnings": [Issue-as-dict...],
        "stats": {...}
    }
    """
    issues = []

    stats = {
        "total_variants": 0,
        "variants_per_contig": {},
        "snp_count": 0,
        "indel_count": 0,
        "multiallelic_count": 0,
        "missing_qual_count": 0,
        "missing_alt_count": 0,
        "invalid_ref_alt_count": 0,
        "nonstandard_bases_count": 0,
        "min_qual": None,
        "max_qual": None,
        "mean_qual": None,
        "qual_n": 0,
        "header_has_contigs": None,
        "sample_count": None,
        "has_required_columns": None,
        "parser_used": None,
        "scanned_records": 0,
        "truncated_scan": False,
    }

    if not os.path.exists(vcf_path):
        issues.append(Issue(
            level="error",
            code="VCF_NOT_FOUND",
            message="VCF file path does not exist.",
            file="vcf",
            where=vcf_path,
            fix="Upload the VCF again or verify the server saved path.",
        ))
        return _finalize(issues, stats)

    # Gzipped VCF requires pysam for parsing
    if vcf_path.lower().endswith(".gz") and not _HAVE_PYSAM:
        issues.append(Issue(
            level="error",
            code="VCF_GZ_REQUIRES_PYSAM",
            message="Gzipped VCF requires pysam for parsing in this prototype.",
            file="vcf",
            fix="Install pysam or upload an uncompressed .vcf file.",
            verify="Decompress .vcf.gz to .vcf and retry."
        ))
        return _finalize(issues, stats)

    # Quick empty-file check
    if os.path.getsize(vcf_path) == 0:
        issues.append(Issue(
            level="error",
            code="VCF_EMPTY",
            message="VCF file is empty.",
            file="vcf",
            fix="Upload a valid VCF containing a header and at least one variant record."
        ))
        return _finalize(issues, stats)

    # Prefer pysam if available
    if _HAVE_PYSAM:
        try:
            stats["parser_used"] = "pysam"
            return _vcf_checks_with_pysam(vcf_path, max_records_scan, issues, stats)
        except Exception as e:
            issues.append(Issue(
                level="error",
                code="VCF_PARSE_FAILED_PYSAM",
                message="VCF could not be parsed with pysam: {}".format(str(e)),
                file="vcf",
                fix="Ensure the file is a valid VCF (tab-separated, correct header). If it is bgzipped (.vcf.gz), ensure it is not corrupted.",
                verify="Try opening the VCF in a viewer or run: bcftools view <file> | head"
            ))
            # fallback to manual scan to provide more actionable hints
            stats["parser_used"] = "manual_fallback"
            return _vcf_checks_manual(vcf_path, max_records_scan, issues, stats)

    # No pysam installed -> manual checks
    stats["parser_used"] = "manual"
    return _vcf_checks_manual(vcf_path, max_records_scan, issues, stats)


def parse_vcf_and_check(vcf_path: str, max_records_scan: int = 200000):
    """
    Returns:
      vcf_info: stats dict
      errors:   list of dicts (title, details, fix, verify, where)
      warnings: list of dicts
      passed:   list of strings
    """
    result = vcf_rule_checks(vcf_path, max_records_scan=max_records_scan)
    stats = result.get("stats", {})
    raw_errors = result.get("errors", [])
    raw_warnings = result.get("warnings", [])

    def _where_obj(where_value):
        where = {"file": "VCF"}
        if where_value:
            text = str(where_value)
            if text.startswith("line "):
                try:
                    where["line"] = int(text.split("line ", 1)[1])
                except Exception:
                    where["line"] = text
            elif ":" in text and text != "NA":
                parts = text.split(":")
                if len(parts) >= 2:
                    where["contig"] = parts[0]
                    try:
                        where["pos"] = int(parts[1])
                    except Exception:
                        where["pos"] = parts[1]
                else:
                    where["line"] = text
            else:
                where["line"] = text
        return where

    def _issue_to_ui(issue, fallback_title):
        title = issue.get("code") or fallback_title
        details = issue.get("message") or ""
        return {
            "title": title,
            "details": details,
            "fix": issue.get("fix"),
            "verify": issue.get("verify"),
            "where": _where_obj(issue.get("where")),
        }

    errors = [_issue_to_ui(e, "VCF error") for e in raw_errors]
    warnings = [_issue_to_ui(w, "VCF warning") for w in raw_warnings]

    passed = []
    if not errors:
        passed.append("VCF file is readable")
    if stats.get("scanned_variants") is not None:
        passed.append("Scanned variants: {}".format(stats.get("scanned_variants", 0)))
    if stats.get("truncated_scan"):
        passed.append("Scan truncated at {} records".format(stats.get("scanned_records", 0)))
    if stats.get("has_required_columns") is True:
        passed.append("Required VCF columns present")

    return stats, errors, warnings, passed


def _vcf_checks_with_pysam(vcf_path: str, max_records_scan: int, issues: List[Issue], stats: Dict[str, Any]) -> Dict[str, Any]:
    v = pysam.VariantFile(vcf_path)

    # Header sanity
    header = v.header
    samples = list(header.samples)
    stats["sample_count"] = len(samples)

    # Required columns exist by VCF spec; but malformed files might still load partially.
    stats["has_required_columns"] = True

    # Header contigs (##contig lines)
    header_contigs = list(header.contigs.keys())
    stats["header_has_contigs"] = (len(header_contigs) > 0)
    if not stats["header_has_contigs"]:
        issues.append(Issue(
            level="warning",
            code="VCF_HEADER_NO_CONTIGS",
            message="VCF header has no ##contig entries. Some tools require contig metadata.",
            file="vcf",
            fix="If needed, regenerate the VCF with a tool that writes contig header lines, or add contig lines based on the reference FASTA index (.fai)."
        ))

    qual_sum = 0.0

    scanned = 0
    for rec in v.fetch():
        scanned += 1
        if scanned > max_records_scan:
            stats["truncated_scan"] = True
            break

        chrom = str(rec.contig) if rec.contig is not None else ""
        pos = int(rec.pos) if rec.pos is not None else -1
        ref = str(rec.ref) if rec.ref is not None else ""
        alts = rec.alts  # tuple or None

        # CHROM
        if not chrom:
            issues.append(Issue(
                level="error",
                code="VCF_EMPTY_CHROM",
                message="Variant record has empty CHROM.",
                file="vcf",
                where="record_index {}".format(scanned)
            ))

        # POS must be >=1
        if pos < 1:
            issues.append(Issue(
                level="error",
                code="VCF_INVALID_POS",
                message="POS must be a positive 1-based integer.",
                file="vcf",
                where="{}:{}".format(chrom, pos if pos != -1 else "NA")
            ))

        # REF basic validation
        if (not ref) or (ref == "."):
            issues.append(Issue(
                level="error",
                code="VCF_INVALID_REF",
                message="REF is missing or '.'.",
                file="vcf",
                where="{}:{}".format(chrom, pos)
            ))
        else:
            if any(c not in VALID_BASES for c in ref.upper()):
                stats["nonstandard_bases_count"] += 1
                issues.append(Issue(
                    level="warning",
                    code="VCF_REF_NONSTANDARD_BASES",
                    message="REF contains non-ACGTN characters (may still be valid depending on pipeline).",
                    file="vcf",
                    where="{}:{}".format(chrom, pos),
                    fix="Prefer A/C/G/T/N. If you expect IUPAC ambiguity codes, confirm your downstream tools accept them."
                ))

        # ALT validation
        if not alts:
            stats["missing_alt_count"] += 1
            issues.append(Issue(
                level="error",
                code="VCF_MISSING_ALT",
                message="ALT is missing.",
                file="vcf",
                where="{}:{}".format(chrom, pos)
            ))
        else:
            # Multiallelic
            if len(alts) > 1:
                stats["multiallelic_count"] += 1

            for alt in alts:
                if alt is None or alt == ".":
                    stats["missing_alt_count"] += 1
                    issues.append(Issue(
                        level="error",
                        code="VCF_INVALID_ALT",
                        message="ALT allele is missing or '.'.",
                        file="vcf",
                        where="{}:{}".format(chrom, pos)
                    ))
                    continue
                # symbolic alleles like <DEL> exist; treat as warning, not invalid
                if alt.startswith("<") and alt.endswith(">"):
                    issues.append(Issue(
                        level="warning",
                        code="VCF_SYMBOLIC_ALT",
                        message="ALT is a symbolic allele (e.g., <DEL>). Some simple checks (like SNP/INDEL) may be inaccurate.",
                        file="vcf",
                        where="{}:{}".format(chrom, pos)
                    ))
                    continue
                # breakend/complex ALT (e.g. N[chr1:123[) - skip SNP/INDEL counting
                if "[" in alt or "]" in alt:
                    issues.append(Issue(
                        level="warning",
                        code="VCF_COMPLEX_ALT",
                        message="ALT appears to be a breakend/complex allele; SNP/INDEL counts may be inaccurate.",
                        file="vcf",
                        where="{}:{}".format(chrom, pos)
                    ))
                    continue

                # sequence alleles
                if any(c not in VALID_BASES for c in alt.upper()):
                    stats["nonstandard_bases_count"] += 1
                    issues.append(Issue(
                        level="warning",
                        code="VCF_ALT_NONSTANDARD_BASES",
                        message="ALT contains non-ACGTN characters (may still be valid depending on pipeline).",
                        file="vcf",
                        where="{}:{}".format(chrom, pos)
                    ))

                # SNP vs INDEL stats (only for non-symbolic ALT)
                vt = _variant_type(ref, alt)
                if vt == "SNP":
                    stats["snp_count"] += 1
                else:
                    stats["indel_count"] += 1

        # QUAL stats
        qual = rec.qual
        if qual is None:
            stats["missing_qual_count"] += 1
        else:
            try:
                q = float(qual)
                stats["qual_n"] += 1
                qual_sum += q
                stats["min_qual"] = q if stats["min_qual"] is None else min(stats["min_qual"], q)
                stats["max_qual"] = q if stats["max_qual"] is None else max(stats["max_qual"], q)
            except Exception:
                stats["missing_qual_count"] += 1
                issues.append(Issue(
                    level="warning",
                    code="VCF_QUAL_NOT_NUMERIC",
                    message="QUAL is present but not numeric.",
                    file="vcf",
                    where="{}:{}".format(chrom, pos)
                ))

        # per-contig counts
        if chrom:
            stats["variants_per_contig"][chrom] = stats["variants_per_contig"].get(chrom, 0) + 1

    stats["scanned_records"] = scanned if scanned <= max_records_scan else max_records_scan
    stats["scanned_variants"] = sum(stats["variants_per_contig"].values())
    stats["total_variants"] = None if stats["truncated_scan"] else stats["scanned_variants"]

    if stats["qual_n"] > 0:
        stats["mean_qual"] = qual_sum / float(stats["qual_n"])

    # High-level checks
    if stats["total_variants"] == 0:
        issues.append(Issue(
            level="warning",
            code="VCF_NO_VARIANTS",
            message="No variant records found (only header).",
            file="vcf",
            fix="Confirm you uploaded the correct file; a VCF should usually contain variant lines after the header."
        ))

    if stats["truncated_scan"]:
        issues.append(Issue(
            level="warning",
            code="VCF_SCAN_TRUNCATED",
            message="Stopped scanning after {} records for performance.".format(max_records_scan),
            file="vcf",
            fix="For large VCFs, this is expected. Rule checks are based on a scan limit."
        ))

    return _finalize(issues, stats)


def _vcf_checks_manual(vcf_path: str, max_records_scan: int, issues: List[Issue], stats: Dict[str, Any]) -> Dict[str, Any]:
    header_found = False
    cols_found = False
    required_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

    qual_sum = 0.0

    scanned = 0
    with open(vcf_path, "r", encoding="utf-8", errors="replace") as f:
        for lineno, line in enumerate(f, start=1):
            line = line.rstrip("\n\r")
            if not line:
                continue

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header_found = True
                parts = line.lstrip("#").split("\t")
                cols_found = all(c in parts for c in required_cols)
                stats["has_required_columns"] = cols_found
                if not cols_found:
                    issues.append(Issue(
                        level="error",
                        code="VCF_MISSING_REQUIRED_COLUMNS",
                        message="VCF header is missing one or more required columns.",
                        file="vcf",
                        where="line {}".format(lineno),
                        fix="Ensure VCF has standard columns: CHROM POS ID REF ALT QUAL FILTER INFO (tab-separated)."
                    ))
                continue

            if not header_found:
                issues.append(Issue(
                    level="error",
                    code="VCF_HEADER_MISSING",
                    message="VCF header line (#CHROM ...) not found before variant records.",
                    file="vcf",
                    where="line {}".format(lineno),
                    fix="Ensure the VCF includes proper header lines beginning with ## and a column header line starting with #CHROM."
                ))
                # keep scanning to give more hints
                header_found = True  # avoid spamming the same error repeatedly

            # Variant line
            scanned += 1
            if scanned > max_records_scan:
                stats["truncated_scan"] = True
                break

            parts = line.split("\t")
            if len(parts) < 8:
                issues.append(Issue(
                    level="error",
                    code="VCF_TOO_FEW_COLUMNS",
                    message="Variant line has fewer than 8 tab-separated columns.",
                    file="vcf",
                    where="line {}".format(lineno),
                    fix="VCF variant lines must have at least 8 columns: CHROM POS ID REF ALT QUAL FILTER INFO."
                ))
                continue

            chrom, pos, _id, ref, alt, qual, _filter, info = parts[:8]

            if chrom == "" or chrom == ".":
                issues.append(Issue("error", "VCF_EMPTY_CHROM", "CHROM is missing.", "vcf", where="line {}".format(lineno)))

            if (not _is_int(pos)) or int(pos) < 1:
                issues.append(Issue("error", "VCF_INVALID_POS", "POS must be a positive integer.", "vcf", where="line {}".format(lineno)))

            if ref == "" or ref == ".":
                issues.append(Issue("error", "VCF_INVALID_REF", "REF is missing or '.'.", "vcf", where="line {}".format(lineno)))
            else:
                if _looks_like_iupac(ref):
                    stats["nonstandard_bases_count"] += 1
                    issues.append(Issue("warning", "VCF_REF_NONSTANDARD_BASES", "REF contains non-ACGTN characters.", "vcf", where="line {}".format(lineno)))

            if alt == "" or alt == ".":
                stats["missing_alt_count"] += 1
                issues.append(Issue("error", "VCF_MISSING_ALT", "ALT is missing or '.'.", "vcf", where="line {}".format(lineno)))
            else:
                alts = alt.split(",")
                if len(alts) > 1:
                    stats["multiallelic_count"] += 1
                for a in alts:
                    a = a.strip()
                    if a.startswith("<") and a.endswith(">"):
                        issues.append(Issue("warning", "VCF_SYMBOLIC_ALT", "ALT is symbolic (e.g., <DEL>).", "vcf", where="line {}".format(lineno)))
                        continue
                    if "[" in a or "]" in a:
                        issues.append(Issue("warning", "VCF_COMPLEX_ALT", "ALT looks like a breakend/complex allele.", "vcf", where="line {}".format(lineno)))
                        continue
                    if a == "" or a == ".":
                        stats["missing_alt_count"] += 1
                        issues.append(Issue("error", "VCF_INVALID_ALT", "ALT allele is missing or '.'.", "vcf", where="line {}".format(lineno)))
                        continue
                    if _looks_like_iupac(a):
                        stats["nonstandard_bases_count"] += 1
                        issues.append(Issue("warning", "VCF_ALT_NONSTANDARD_BASES", "ALT contains non-ACGTN characters.", "vcf", where="line {}".format(lineno)))

                    vt = _variant_type(ref, a)
                    if vt == "SNP":
                        stats["snp_count"] += 1
                    else:
                        stats["indel_count"] += 1

            if qual == "." or qual == "":
                stats["missing_qual_count"] += 1
            else:
                if _is_float(qual):
                    q = float(qual)
                    stats["qual_n"] += 1
                    qual_sum += q
                    stats["min_qual"] = q if stats["min_qual"] is None else min(stats["min_qual"], q)
                    stats["max_qual"] = q if stats["max_qual"] is None else max(stats["max_qual"], q)
                else:
                    stats["missing_qual_count"] += 1
                    issues.append(Issue("warning", "VCF_QUAL_NOT_NUMERIC", "QUAL is present but not numeric.", "vcf", where="line {}".format(lineno)))

            # per-contig counts
            if chrom:
                stats["variants_per_contig"][chrom] = stats["variants_per_contig"].get(chrom, 0) + 1

    stats["scanned_records"] = scanned if scanned <= max_records_scan else max_records_scan
    stats["scanned_variants"] = sum(stats["variants_per_contig"].values())
    stats["total_variants"] = None if stats["truncated_scan"] else stats["scanned_variants"]
    if stats["qual_n"] > 0:
        stats["mean_qual"] = qual_sum / float(stats["qual_n"])

    if not header_found:
        issues.append(Issue(
            level="error",
            code="VCF_HEADER_MISSING",
            message="VCF header not found.",
            file="vcf",
            fix="Ensure the file contains '##' meta lines and '#CHROM' header line."
        ))

    if stats["truncated_scan"]:
        issues.append(Issue(
            level="warning",
            code="VCF_SCAN_TRUNCATED",
            message="Stopped scanning after {} records for performance.".format(max_records_scan),
            file="vcf"
        ))

    return _finalize(issues, stats)


def _finalize(issues: List[Issue], stats: Dict[str, Any]) -> Dict[str, Any]:
    errors = [i.to_dict() for i in issues if i.level == "error"]
    warnings = [i.to_dict() for i in issues if i.level == "warning"]
    return {"errors": errors, "warnings": warnings, "stats": stats}
