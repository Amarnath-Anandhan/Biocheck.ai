from typing import Dict, List, Tuple

ALLOWED_DNA = set("ACGTN")  # prototype: only DNA + N


def parse_fasta_and_check(fasta_path: str) -> Tuple[Dict[str, int], List[dict], List[dict], List[str]]:
    """
    Returns:
      contig_lengths: {contig_name: length}
      errors:   list of dicts (title, details, fix, verify, where)
      warnings: list of dicts
      passed:   list of strings
    """

    errors: List[dict] = []
    warnings: List[dict] = []
    passed: List[str] = []

    contig_lengths: Dict[str, int] = {}

    current_name = None
    current_len = 0
    seen_names = set()

    def finalize_record(line_no: int):
        nonlocal current_name, current_len
        if current_name is None:
            return
        if current_len == 0:
            errors.append({
                "title": "Empty sequence",
                "details": f"Contig '{current_name}' has no sequence lines.",
                "fix": "Remove the empty contig or add its sequence.",
                "verify": "Confirm there are DNA letters after the contig header.",
                "where": {"file": "FASTA", "contig": current_name, "line": line_no}
            })
        else:
            contig_lengths[current_name] = current_len

        current_name = None
        current_len = 0

    try:
        f = open(fasta_path, "r", encoding="utf-8", errors="replace")
    except Exception as e:
        errors.append({
            "title": "FASTA not readable",
            "details": "Could not read file: {}".format(e),
            "fix": "Upload a valid FASTA text file.",
            "verify": "Open the FASTA in a text editor and confirm it contains '>' headers and sequences.",
            "where": {"file": "FASTA"}
        })
        return {}, errors, warnings, passed

    seen_any = False
    first_nonempty = None
    line_no = 0
    for raw in f:
        line_no += 1
        line = raw.rstrip("\n")

        if not line.strip():
            continue

        if not seen_any:
            seen_any = True
            first_nonempty = (line_no, line)
            if not line.startswith(">"):
                errors.append({
                    "title": "Invalid FASTA start",
                    "details": "FASTA must begin with a header line starting with '>' (first non-empty line is not a header).",
                    "fix": "Ensure the first sequence has a header like >chr1.",
                    "verify": "Open file and confirm the first non-empty line starts with '>'.",
                    "where": {"file": "FASTA", "line": line_no}
                })

        if line.startswith(">"):
            # finalize previous record
            finalize_record(line_no)

            header = line[1:].strip()
            if header == "":
                errors.append({
                    "title": "Missing contig name",
                    "details": "A header line '>' exists but contig name is empty.",
                    "fix": "Give a contig name like >chr1.",
                    "verify": "Confirm each header line has a name after '>'.",
                    "where": {"file": "FASTA", "line": line_no}
                })
                current_name = None
                continue

            name = header.split()[0]  # take first token as contig id

            if name in seen_names:
                errors.append({
                    "title": "Duplicate contig name",
                    "details": f"Contig name '{name}' appears more than once.",
                    "fix": "Rename one of the duplicate headers or remove duplicates.",
                    "verify": "Confirm each header name is unique.",
                    "where": {"file": "FASTA", "contig": name, "line": line_no}
                })
            else:
                seen_names.add(name)

            current_name = name
            current_len = 0

        else:
            if current_name is None:
                errors.append({
                    "title": "Sequence without header",
                    "details": "Found sequence line before any valid header.",
                    "fix": "Add a header line starting with '>' before the sequence.",
                    "verify": "Confirm all sequences are under a header.",
                    "where": {"file": "FASTA", "line": line_no}
                })
                continue

            seq = line.strip().upper().replace(" ", "").replace("\t", "")
            # check invalid characters (warn, since some pipelines allow IUPAC)
            bad = {ch for ch in seq if ch not in ALLOWED_DNA}
            if bad:
                warnings.append({
                    "title": "Non-standard bases in sequence",
                    "details": "Contig '{}' contains non-ACGTN characters: {}".format(current_name, ", ".join(sorted(bad))),
                    "fix": "If this is a DNA FASTA, consider replacing with A/C/G/T/N.",
                    "verify": "Confirm downstream tools accept IUPAC ambiguity codes.",
                    "where": {"file": "FASTA", "contig": current_name, "line": line_no}
                })

            current_len += len(seq)

    f.close()

    if not seen_any:
        errors.append({
            "title": "Empty FASTA file",
            "details": "FASTA file has no content.",
            "fix": "Upload a non-empty FASTA file.",
            "verify": "Confirm FASTA has at least one header line starting with '>'.",
            "where": {"file": "FASTA"}
        })
        return {}, errors, warnings, passed

    # finalize last record
    finalize_record(line_no=line_no)

    # Passed checks (only if not major failures)
    if not errors:
        passed.append("FASTA file is readable")
        passed.append("FASTA headers and sequences parsed successfully")
        passed.append(f"Total contigs: {len(contig_lengths)}")

    # Warnings (optional, light)
    if contig_lengths and len(contig_lengths) > 2000:
        warnings.append({
            "title": "Large FASTA detected",
            "details": f"FASTA contains {len(contig_lengths)} contigs; processing may be slower in prototype.",
            "fix": "For testing, use a smaller FASTA (one chromosome or small genome).",
            "verify": "Try again with a smaller FASTA and confirm results load quickly.",
            "where": {"file": "FASTA"}
        })

    return contig_lengths, errors, warnings, passed
