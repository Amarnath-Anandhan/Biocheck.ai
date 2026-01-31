# Manual test runner for BioCheckAI
import os
import tempfile
from rules.fasta_checks import parse_fasta_and_check
from rules.vcf_checks import parse_vcf_and_check
from rules.cross_checks import cross_check_fasta_vcf


def _write(tmpdir, name, content):
    path = os.path.join(tmpdir, name)
    with open(path, 'w', encoding='utf-8') as f:
        f.write(content)
    return path


def run():
    failures = 0
    with tempfile.TemporaryDirectory() as tmp:
        # FASTA tests
        fasta_valid = ">chr1\nACGTACGT\n>chr2\nNNNN\n"
        fasta_empty = ""
        fasta_dup = ">chr1\nACGT\n>chr1\nACGT\n"
        fasta_invalid = ">chr1\nACGTXYZ\n"

        tests = [
            ("FASTA valid", fasta_valid, False),
            ("FASTA empty", fasta_empty, True),
            ("FASTA duplicate contig", fasta_dup, True),
            ("FASTA invalid chars", fasta_invalid, True),
        ]

        for name, content, expect_error in tests:
            path = _write(tmp, "test.fasta", content)
            contigs, errors, warnings, passed = parse_fasta_and_check(path)
            ok = (len(errors) == 0)
            if ok == (not expect_error):
                print("PASS:", name)
            else:
                print("FAIL:", name, "errors=", len(errors))
                failures += 1

        # VCF tests
        vcf_valid = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1\t.\tA\tT\t50\tPASS\t.
"""
        vcf_no_chrom = """##fileformat=VCFv4.2
chr1\t1\t.\tA\tT\t50\tPASS\t.
"""
        vcf_space_sep = """##fileformat=VCFv4.2
#CHROM POS ID REF ALT QUAL FILTER INFO
chr1 1 . A T 50 PASS .
"""
        vcf_pos0 = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t0\t.\tA\tT\t50\tPASS\t.
"""
        vcf_alt_dot = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1\t.\tA\t.\t50\tPASS\t.
"""
        vcf_multialt = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t2\t.\tA\tT,G\t50\tPASS\t.
"""
        vcf_symbolic = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t3\t.\tA\t<DEL>\t50\tPASS\t.
"""

        vcf_tests = [
            ("VCF valid", vcf_valid, False),
            ("VCF missing #CHROM", vcf_no_chrom, True),
            ("VCF space separated", vcf_space_sep, True),
            ("VCF POS=0", vcf_pos0, True),
            ("VCF ALT=.", vcf_alt_dot, True),
            ("VCF multiallelic", vcf_multialt, False),
            ("VCF symbolic ALT", vcf_symbolic, False),
        ]

        for name, content, expect_error in vcf_tests:
            path = _write(tmp, "test.vcf", content)
            stats, errors, warnings, passed = parse_vcf_and_check(path)
            ok = (len(errors) == 0)
            if ok == (not expect_error):
                print("PASS:", name)
            else:
                print("FAIL:", name, "errors=", len(errors))
                failures += 1

        # Cross-check tests
        fasta_path = _write(tmp, "cross.fasta", fasta_valid)
        contigs, errors, warnings, passed = parse_fasta_and_check(fasta_path)

        vcf_missing_contig = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chrX\t1\t.\tA\tT\t50\tPASS\t.
"""
        vcf_out_of_range = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t999\t.\tA\tT\t50\tPASS\t.
"""
        vcf_chr_mismatch = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t1\t.\tA\tT\t50\tPASS\t.
"""

        for name, content in [
            ("Cross missing contig", vcf_missing_contig),
            ("Cross out of range", vcf_out_of_range),
            ("Cross chr mismatch", vcf_chr_mismatch),
        ]:
            vcf_path = _write(tmp, "cross.vcf", content)
            cross = cross_check_fasta_vcf(contigs, vcf_path, ref_check=False)
            if cross.get("errors"):
                print("PASS:", name)
            else:
                print("FAIL:", name)
                failures += 1

    return failures


if __name__ == "__main__":
    fails = run()
    raise SystemExit(1 if fails else 0)
