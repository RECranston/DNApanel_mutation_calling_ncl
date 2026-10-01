#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --partition=default_free
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --job-name=script6_merge_phased
#SBATCH --output=logs/merge_phased_%A_%a.out
#SBATCH --array=1-100%100

# Script to merge phased Mutect2 VCF records (sharing a PID) into a single minimal
# compound variant record before VEP annotation. VEP annotates each input record
# independently and does not use phasing to combine nearby variants into one
# consequence prediction (confirmed limitation, see Ensembl VEP paper, Genome Biology
# 2016) - so several small linked indels/substitutions that together form one real
# in-frame event can each look like a spurious frameshift if left unmerged.


# Ruth Cranston 2026
[ $# -ne 3 ] && { echo -en \
"\nRuth Cranston 2026\n\n
*** Script to merge phased Mutect2 records into minimal compound variants, ready for VEP.
Runs in current directory. Input dir is location of PASS-filtered Mutect2 VCFs (script5 output).
Output directory is created.
<sample sheet> <input dir (relative)> <output dir (relative)>
example run: sbatch ./script6_merge_phased_variants.sh sample_sheet.txt output_mutation_calling/ output_merged_phased/ *** \n\n" ; exit 1; }

# Set variables
BASE_DIR="$PWD"
ASSEMBLY="GRCh37"
SAMPLE_SHEET=$1
INPUT_DIR=${BASE_DIR}/$2
OUTPUT_DIR=${BASE_DIR}/$3
REFERENCE_DIR=${BASE_DIR}/References/${ASSEMBLY}

# Max bp span (first to last record) allowed within a PID group before merging - matches
# GATK's own default "bad haplotype" phasing distance. Lower this (e.g. to 10-20) if you'd
# rather only merge tightly-clustered indels and leave more distant phased substitutions
# annotated separately.
MAX_PHASED_SPAN=10

echo -en " * Loading modules...\n"
module --force purge
module load Python/3.11.5-GCCcore-13.2.0
module load SAMtools/1.21-GCC-13.3.0
module load HTSlib/1.21-GCC-13.3.0
set -euo pipefail
echo -en " * Environment set up.\n"

mkdir -p ${OUTPUT_DIR}
mkdir -p logs

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET})
SAMPLE_ID=$(echo $LINE | awk '{print $1}')
if [[ -z "${SAMPLE_ID}" ]]; then
    echo "ERROR: SAMPLE_ID is empty for array task ${SLURM_ARRAY_TASK_ID} - check row ${SLURM_ARRAY_TASK_ID} of ${SAMPLE_SHEET} exists and is correctly formatted" >&2
    exit 1
fi
echo "Processing sample: ${SAMPLE_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"

if [ "${ASSEMBLY}" == "GRCh38" ]; then
    REF_FASTA=${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta
else
    REF_FASTA=${REFERENCE_DIR}/Homo_sapiens_assembly19.fasta
fi

INPUT_VCF=${INPUT_DIR}${SAMPLE_ID}_tumor_filtered_PASS.vcf.gz
MERGED_VCF=${OUTPUT_DIR}${SAMPLE_ID}_tumor_filtered_PASS_merged.vcf.gz

if [[ ! -f "${INPUT_VCF}" ]]; then
    echo "ERROR: expected input VCF not found at ${INPUT_VCF}" >&2
    exit 1
fi

MERGE_SCRIPT="${OUTPUT_DIR}.merge_phased_${SLURM_ARRAY_TASK_ID}.py"
cat > ${MERGE_SCRIPT} << 'PYEOF'
#!/usr/bin/env python3
"""
Merge phased VCF records sharing a PID into a single minimal delins record,
so downstream annotation (VEP) sees one clean compound variant instead of
several independent small edits that can each look like a frameshift.

Ruth Cranston 2026
"""
import sys
import gzip
import subprocess
import argparse


def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def fetch_ref_window(ref_fasta, chrom, start, end):
    """Fetch reference bases [start, end] (1-based, inclusive) via samtools faidx."""
    region = f"{chrom}:{start}-{end}"
    out = subprocess.run(
        ["samtools", "faidx", ref_fasta, region],
        capture_output=True, text=True, check=True
    ).stdout
    lines = out.strip().splitlines()
    return "".join(lines[1:])  # drop the ">chr:start-end" header line


def get_pid(fmt_field, sample_field):
    fmt_keys = fmt_field.split(":")
    sample_vals = sample_field.split(":")
    if "PID" not in fmt_keys:
        return None
    idx = fmt_keys.index("PID")
    if idx >= len(sample_vals):
        return None
    val = sample_vals[idx]
    return val if val not in (".", "") else None


def reconstruct_group(ref_fasta, chrom, group):
    """group: list of (pos, ref, alt) tuples, sorted by pos."""
    group_start = group[0][0]
    group_end = max(pos + len(ref) - 1 for pos, ref, alt in group)

    ref_span = fetch_ref_window(ref_fasta, chrom, group_start, group_end)
    ref_bases = {group_start + i: b for i, b in enumerate(ref_span)}
    variants = {pos: (ref, alt) for pos, ref, alt in group}

    i = group_start
    seg = []
    while i <= group_end:
        if i in variants:
            ref, alt = variants[i]
            seg.append(alt)
            i += len(ref)
        else:
            seg.append(ref_bases[i])
            i += 1
    mutant_span = "".join(seg)

    # Trim to the minimal representation, but NEVER let REF or ALT become empty - VCF
    # requires at least one base in each. Trim common suffix first, then common prefix,
    # sequentially (not independently), each time stopping once either string would drop
    # to length 1. This is standard VCF left-normalization with an anchor base preserved.
    ref, alt = ref_span, mutant_span
    minimal_pos = group_start

    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]

    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        minimal_pos += 1

    if not ref or not alt:
        raise ValueError(f"trimming produced an empty REF or ALT (ref={ref!r}, alt={alt!r})")

    return minimal_pos, ref, alt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-vcf", required=True)
    ap.add_argument("--output-vcf", required=True)
    ap.add_argument("--ref-fasta", required=True)
    ap.add_argument("--max-span", type=int, default=100)
    args = ap.parse_args()

    header_lines = []
    data_lines = []
    with open_maybe_gzip(args.input_vcf) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header_lines.append(line)
            elif line:
                data_lines.append(line.split("\t"))

    groups = []
    current_key = None
    current_group = []
    for rec in data_lines:
        chrom, pos, vid, ref, alt, qual, filt, info = rec[:8]
        rest = rec[8:]
        pid = get_pid(rest[0], rest[1]) if len(rest) >= 2 else None
        key = (chrom, pid)
        if pid is not None and key == current_key:
            current_group.append(rec)
        else:
            if current_group:
                groups.append((current_key[0], current_key[1], current_group))
            current_group = [rec]
            current_key = key
    if current_group:
        groups.append((current_key[0], current_key[1], current_group))

    merged_count = 0
    out_lines = []
    for chrom, pid, recs in groups:
        if pid is None or len(recs) < 2:
            out_lines.extend(recs)
            continue

        # Within a PID group, only merge records that are strictly adjacent - i.e. the next
        # record's POS is exactly where the previous record's REF ends, with zero bases of
        # untouched reference in between. This distinguishes a genuine compound indel (several
        # small edits abutting each other, which Mutect2 sometimes splits across records) from
        # a separately-phased but more distant substitution that happens to share a PID -
        # sharing a PID only means "same allele", not "the same compound event".
        subclusters = [[recs[0]]]
        for rec in recs[1:]:
            prev = subclusters[-1][-1]
            prev_pos, prev_ref = int(prev[1]), prev[3]
            cur_pos = int(rec[1])
            if cur_pos == prev_pos + len(prev_ref):
                subclusters[-1].append(rec)
            else:
                subclusters.append([rec])

        for sub in subclusters:
            if len(sub) < 2:
                out_lines.extend(sub)
                continue

            positions = [int(r[1]) for r in sub]
            span = max(positions) - min(positions)
            if span > args.max_span:
                out_lines.extend(sub)
                continue

            group = [(int(r[1]), r[3], r[4]) for r in sub]
            try:
                new_pos, new_ref, new_alt = reconstruct_group(args.ref_fasta, chrom, group)
            except Exception as e:
                sys.stderr.write(f"WARNING: could not merge group at {chrom}:{positions} (PID={pid}): {e}\n")
                out_lines.extend(sub)
                continue

            base = sub[0]
            merged = list(base)
            merged[1] = str(new_pos)
            merged[3] = new_ref
            merged[4] = new_alt
            merged[7] = base[7] + f";PHASED_MERGE_COUNT={len(sub)}"
            out_lines.append(merged)
            merged_count += 1
            sys.stderr.write(f"Merged {len(sub)} records at {chrom}:{positions} (PID={pid}) -> {chrom}:{new_pos} {new_ref}>{new_alt}\n")

    new_info_header = '##INFO=<ID=PHASED_MERGE_COUNT,Number=1,Type=Integer,Description="Number of phased records merged into this compound variant">'
    with open_maybe_gzip(args.output_vcf, "wt") as out:
        for h in header_lines:
            if h.startswith("#CHROM"):
                out.write(new_info_header + "\n")
            out.write(h + "\n")
        for rec in out_lines:
            out.write("\t".join(rec) + "\n")

    sys.stderr.write(f"Done. Merged {merged_count} phased group(s).\n")


if __name__ == "__main__":
    main()
PYEOF

python3 ${MERGE_SCRIPT} \
    --input-vcf ${INPUT_VCF} \
    --output-vcf ${MERGED_VCF%.gz} \
    --ref-fasta ${REF_FASTA} \
    --max-span ${MAX_PHASED_SPAN}
echo -ne "*** Phased-variant merge complete! ***\n"

bgzip ${MERGED_VCF%.gz}
tabix -p vcf ${MERGED_VCF}
rm -f ${MERGE_SCRIPT}
echo -ne "*** Indexing complete! ***\n"
echo -ne "*** All done! ***\n"
