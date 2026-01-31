# BioCheckAI (Prototype)

BioCheckAI is a Flask-based prototype for validating FASTA and VCF files and reporting consistency issues, errors, and warnings in a user-friendly report.

## Features
- Upload FASTA/VCF files or paste content directly.
- FASTA validation (headers, duplicates, invalid characters).
- VCF validation (required columns, POS, REF/ALT, QUAL).
- Cross-checks between FASTA and VCF (contig existence, out-of-range POS, naming mismatch).
- Grouped error reporting for large files.

## Requirements
- Python 3.12 recommended (older versions may work but are not guaranteed).
- Flask
- Optional: `pysam` for better VCF/FASTA parsing and gzipped VCF support.

## Setup
```powershell
python -m pip install flask
```

Optional (for gzipped VCF / indexed FASTA support):
```powershell
python -m pip install pysam
```

## Run
```powershell
python app.py
```
Then open:
```
http://127.0.0.1:5000
```

## Manual Tests
```powershell
python tests/manual_test_runner.py
```

## Notes / Limitations
- Large uploads are capped by `MAX_CONTENT_LENGTH` in `app.py`.
- Gzipped VCF requires `pysam`. Without it, `.vcf.gz` files are rejected.
- SNP/INDEL counts are simplified and may not reflect complex/breakend alleles.

## License
This repository is proprietary and **all rights are reserved**. See `LICENSE`.
