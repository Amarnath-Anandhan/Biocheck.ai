import os
import uuid
from flask import Flask, render_template, request
from werkzeug.utils import secure_filename
from rules.cross_checks import cross_check_fasta_vcf
from rules.fasta_checks import parse_fasta_and_check
from rules.vcf_checks import parse_vcf_and_check
from rules.cross_checks import cross_check_fasta_vcf

app = Flask(__name__)

# Folder to store uploaded files
UPLOAD_FOLDER = os.path.join(app.root_path, "uploads")
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER
app.config["MAX_CONTENT_LENGTH"] = 50 * 1024 * 1024  # 50 MB upload limit


@app.errorhandler(413)
def request_entity_too_large(_error):
    report = {
        "summary": {"status": "ERROR", "error_count": 1, "warning_count": 0},
        "errors": [{
            "title": "Upload too large",
            "details": "The uploaded file exceeds the server limit.",
            "fix": "Use a smaller file or increase MAX_CONTENT_LENGTH in app.py.",
            "verify": "Try again with a smaller file size.",
            "where": {"file": "UPLOAD"},
        }],
        "warnings": [],
        "passed": [],
    }
    return render_template("result.html", report=report), 413


@app.get("/")
def home():
    return render_template("home.html")


@app.get("/upload")
def upload_page():
    return render_template("uploads.html")


@app.post("/check")
def run_check():
    # 1) Receive files
    fasta_file = request.files.get("fasta_file")
    vcf_file = request.files.get("vcf_file")
    fasta_text = (request.form.get("fasta_text") or "").strip()
    vcf_text = (request.form.get("vcf_text") or "").strip()

    # 2) Basic validations (prototype-level)
    errors = []
    has_fasta_file = fasta_file is not None and fasta_file.filename.strip() != ""
    has_vcf_file = vcf_file is not None and vcf_file.filename.strip() != ""

    if not has_fasta_file and not fasta_text:
        errors.append({
            "title": "Missing FASTA file",
            "details": "Please upload a FASTA file.",
            "fix": "Choose a .fasta/.fa file and try again.",
            "verify": "Confirm the FASTA file is selected before clicking Run Check."
        })

    if not has_vcf_file and not vcf_text:
        errors.append({
            "title": "Missing VCF file",
            "details": "Please upload a VCF file.",
            "fix": "Choose a .vcf file and try again.",
            "verify": "Confirm the VCF file is selected before clicking Run Check."
        })

    # If missing files, show result immediately
    if errors:
        report = {
            "summary": {"status": "ERROR", "error_count": len(errors), "warning_count": 0},
            "errors": errors,
            "warnings": [],
            "passed": []
        }
        return render_template("result.html", report=report)

    # 3) Save files or pasted text
    if has_fasta_file:
        fasta_name = secure_filename(fasta_file.filename) or "upload.fasta"
        fasta_name = "{}_{}".format(uuid.uuid4().hex, fasta_name)
        fasta_path = os.path.join(app.config["UPLOAD_FOLDER"], fasta_name)
        fasta_file.save(fasta_path)
    else:
        fasta_name = "{}_pasted_fasta.fasta".format(uuid.uuid4().hex)
        fasta_path = os.path.join(app.config["UPLOAD_FOLDER"], fasta_name)
        with open(fasta_path, "w", encoding="utf-8") as handle:
            handle.write(fasta_text + "\n")

    if has_vcf_file:
        vcf_name = secure_filename(vcf_file.filename) or "upload.vcf"
        vcf_name = "{}_{}".format(uuid.uuid4().hex, vcf_name)
        vcf_path = os.path.join(app.config["UPLOAD_FOLDER"], vcf_name)
        vcf_file.save(vcf_path)
    else:
        vcf_name = "{}_pasted_vcf.vcf".format(uuid.uuid4().hex)
        vcf_path = os.path.join(app.config["UPLOAD_FOLDER"], vcf_name)
        with open(vcf_path, "w", encoding="utf-8") as handle:
            handle.write(vcf_text + "\n")

    # 4) FASTA + VCF checks
    contig_lengths, fasta_errors, fasta_warnings, fasta_passed = parse_fasta_and_check(fasta_path)
    vcf_info, vcf_errors, vcf_warnings, vcf_passed = parse_vcf_and_check(vcf_path)

    errors = fasta_errors + vcf_errors
    warnings = fasta_warnings + vcf_warnings
    passed = fasta_passed + vcf_passed

    # 5) Cross-checks (contig existence / range / naming mismatch)
    cross = cross_check_fasta_vcf(contig_lengths, vcf_path, ref_check=False)

    def _normalize_issue(issue):
        if not isinstance(issue, dict):
            return issue
        if "title" in issue and "details" in issue:
            return issue
        where = {"file": "CROSS"}
        if issue.get("where"):
            where["line"] = issue.get("where")
        return {
            "title": issue.get("code") or "CROSS_CHECK",
            "details": issue.get("message") or "",
            "fix": issue.get("fix"),
            "verify": issue.get("verify"),
            "where": where,
        }

    errors += [_normalize_issue(e) for e in cross.get("errors", [])]
    warnings += [_normalize_issue(w) for w in cross.get("warnings", [])]

    def _cap_issues(items, max_items, label):
        if len(items) <= max_items:
            return items
        remaining = len(items) - max_items
        capped = items[:max_items]
        capped.append({
            "title": "{} more {} not shown".format(remaining, label),
            "details": "Too many similar records. Showing first {}.".format(max_items),
            "where": {"file": "SYSTEM"},
        })
        return capped

    # Cap error/warning lists to keep UI responsive
    errors = _cap_issues(errors, 200, "errors")
    warnings = _cap_issues(warnings, 200, "warnings")

    def _compress_numbers(nums, max_ranges=8):
        nums = sorted({n for n in nums if isinstance(n, int)})
        ranges = []
        if not nums:
            return ranges
        start = prev = nums[0]
        for n in nums[1:]:
            if n == prev + 1:
                prev = n
                continue
            ranges.append((start, prev))
            start = prev = n
        ranges.append((start, prev))
        if len(ranges) > max_ranges:
            ranges = ranges[:max_ranges] + [("...", "...")]
        return ranges

    def _group_issues(items):
        grouped = {}
        for item in items:
            key = (
                item.get("title"),
                item.get("details"),
                item.get("fix"),
                item.get("verify"),
                item.get("where", {}).get("file"),
            )
            entry = grouped.get(key)
            if entry is None:
                entry = {
                    "title": item.get("title"),
                    "details": item.get("details"),
                    "fix": item.get("fix"),
                    "verify": item.get("verify"),
                    "where": {
                        "file": item.get("where", {}).get("file"),
                        "lines": [],
                    },
                    "_count": 0,
                }
                grouped[key] = entry
            entry["_count"] += 1
            where = item.get("where") or {}
            if "line" in where:
                entry["where"]["lines"].append(where.get("line"))
            elif "pos" in where and "contig" in where:
                # keep example positions only
                entry.setdefault("examples", [])
                if len(entry["examples"]) < 5:
                    entry["examples"].append("{}:{}".format(where.get("contig"), where.get("pos")))

        result = []
        for entry in grouped.values():
            lines = entry["where"].get("lines") or []
            entry["where"]["line_ranges"] = _compress_numbers(lines)
            if entry["_count"] > 1:
                entry["details"] = "{} ({} occurrences)".format(entry["details"], entry["_count"])
            entry.pop("_count", None)
            result.append(entry)
        return result

    errors = _group_issues(errors)
    warnings = _group_issues(warnings)

    status = "OK"
    if errors:
        status = "ERROR"
    elif warnings:
        status = "WARNING"

    report = {
        "summary": {
            "status": status,
            "error_count": len(errors),
            "warning_count": len(warnings),
        },
        "errors": errors,
        "warnings": warnings,
        "passed": passed,
        "cross": cross,
        "debug": {
            "fasta_contigs": len(contig_lengths),
            "vcf_scanned_variants": vcf_info.get("scanned_variants", 0),
            "fasta_path": fasta_path,
            "vcf_path": vcf_path,
        },
    }

    return render_template("result.html", report=report)


if __name__ == "__main__":
    # debug=True for development only
    app.run(debug=True)
