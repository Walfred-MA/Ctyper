#!/usr/bin/env python3

import argparse
import gzip
import re
import sys
from pathlib import Path
import collections as cl

NOTE_PATTERN = re.compile(r"\(([^)]+)\)")


def warning(message):
    """Print a warning without contaminating the output TSV."""
    print(f"WARNING: {message}", file=sys.stderr)


def open_text(path):
    """Open plain-text or gzipped files transparently."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def normalize_build(build):
    """Normalize common genome build names so hg38 == GRCh38 and hg19 == GRCh37."""
    if build is None:
        return None

    cleaned = str(build).strip()
    lower = cleaned.lower()

    if lower in {"hg38", "grch38", "grch38.p13", "grch38.p14"} or lower.startswith("grch38"):
        return "GRCh38"
    if lower in {"hg19", "grch37", "grch37.p13"} or lower.startswith("grch37"):
        return "GRCh37"

    return cleaned


def detect_gencode_metadata(gff3_path):
    """
    Try to detect GENCODE release/version and genome build from the GFF3 header
    and, as a fallback, from the annotation filename.
    """
    metadata = {
        "version": None,
        "genome_build": None,
    }

    header_lines = []
    try:
        with open_text(gff3_path) as f:
            for line in f:
                if line.startswith("#"):
                    header_lines.append(line.strip())
                    continue
                break
    except OSError as exc:
        warning(f"could not read GENCODE annotation header from '{gff3_path}': {exc}")
        return metadata

    header_text = "\n".join(header_lines)
    file_name = Path(gff3_path).name

    version_patterns = [
        r"GENCODE[^0-9vV]*(?:v|version\s*)?([0-9]+)",
        r"gencode\.v([0-9]+)",
        r"\bversion\s+([0-9]+)\b",
    ]
    for pattern in version_patterns:
        match = re.search(pattern, header_text, flags=re.IGNORECASE)
        if match:
            metadata["version"] = match.group(1)
            break

    if metadata["version"] is None:
        match = re.search(r"gencode\.v([0-9]+)", file_name, flags=re.IGNORECASE)
        if match:
            metadata["version"] = match.group(1)

    build_patterns = [
        r"\b(GRCh3[78](?:\.p[0-9]+)?)\b",
        r"\b(hg(?:19|38))\b",
    ]
    for pattern in build_patterns:
        match = re.search(pattern, header_text, flags=re.IGNORECASE)
        if match:
            metadata["genome_build"] = normalize_build(match.group(1))
            break

    if metadata["genome_build"] is None:
        for pattern in build_patterns:
            match = re.search(pattern, file_name, flags=re.IGNORECASE)
            if match:
                metadata["genome_build"] = normalize_build(match.group(1))
                break

    return metadata


def gencode_check_message(expected_version):
    """Message requested for missing-ID warnings."""
    return f"checking your genecode version, default version: {expected_version}"


def validate_gencode_metadata(args, metadata):
    """
    Warn if the detected annotation version/build does not match expectations.
    By default this is warning-only. It exits only when --strict-gencode-version is used.
    """
    expected_version = str(args.expected_gencode_version) if args.expected_gencode_version is not None else None
    detected_version = metadata.get("version")
    expected_build = normalize_build(args.expected_genome_build)
    detected_build = normalize_build(metadata.get("genome_build"))

    mismatches = []


    if expected_build and detected_build and detected_build != expected_build:
        mismatches.append(f"genome build mismatch: expected {args.expected_genome_build}, detected {metadata.get('genome_build')}")
    elif expected_build and not detected_build:
        mismatches.append(f"could not detect genome build from annotation; expected default build {args.expected_genome_build}")

    for message in mismatches:
        warning(f"{message}; {gencode_check_message(expected_version)}")

    if args.strict_gencode_version and mismatches:
        raise SystemExit("Stopping because --strict-gencode-version was used and the annotation did not match expectations.")


def interval_intersection(X, Y, frac_threshold=0.5):

    i = j = 0
    result = []
    x_significant = set()
    X = sorted(X)
    Y = sorted(Y)

    while i < len(X) and j < len(Y):
        a_start, a_end = X[i][:2]
        b_start, b_end = Y[j][:2]

        # Compute overlap
        start = max(a_start, b_start)
        end = min(a_end, b_end)

        if start < end:  # Non-zero intersection
            result.append([start, end])

            # Check fraction of X[i] that is covered
            x_len = a_end - a_start
            overlap_len = end - start
            if x_len > 0 and (overlap_len / x_len) >= frac_threshold:
                x_significant.add(i)

        # Advance the interval that ends first
        if a_end < b_end:
            i += 1
        else:
            j += 1

    return result, x_significant


def compress_ranges(numbers):

    if len(numbers) == 0:
        return ""

    compressed = []
    start = prev = numbers[0]

    for n in numbers[1:]:
        if n == prev + 1:
            prev = n
        else:
            if start == prev:
                compressed.append(f"{start}")
            else:
                compressed.append(f"{start}-{prev}")
            start = prev = n

    # Add final group
    if start == prev:
        compressed.append(f"{start}")
    else:
        compressed.append(f"{start}-{prev}")

    return compressed


def parse_attributes(attr_str):
    """Parse GFF3 attributes column into a dict."""
    attrs = {}
    for entry in attr_str.strip().split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            attrs[key] = value
    return attrs


def readgff3(inf, records):

    transsize = cl.defaultdict(int)
    genetoanno = cl.defaultdict(list)
    transtogenes = cl.defaultdict(str)
    transtoexons = cl.defaultdict(list)
    manelist = cl.defaultdict(lambda: [[], []])
    exon_ids = set()

    with open_text(inf) as f:

        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) != 9 or cols[1] != "HAVANA":
                continue

            chrom, source, feature, start, end, score, strand, phase, attributes = cols
            start, end = int(start), int(end)

            attrs = parse_attributes(attributes)

            if feature == "gene" and "ID" in attrs and "gene_type" in attrs:
                gene_id = attrs["gene_id"].split(".")[0]
                anno = attrs["gene_type"]
                name = attrs["gene_name"] if "gene_name" in attrs else ""

                genetoanno[gene_id] = [name, anno]

            if feature == "transcript":
                gene_id = attrs["gene_id"].split(".")[0]
                enst_id = attrs["transcript_id"].split(".")[0]

                if "MANE_Select" in attrs:
                    manelist[gene_id][0].append(enst_id)
                else:
                    manelist[gene_id][1].append(enst_id)

            if feature == "exon" and "ID" in attrs and "transcript_name" in attrs and "gene_id" in attrs:
                ense_id = attrs["exon_id"].split(".")[0]
                enst_id = attrs["transcript_id"].split(".")[0]
                ensg_id = attrs["gene_id"].split(".")[0]

                exon_ids.add(ense_id)
                transtogenes[enst_id] = ensg_id
                transtoexons[enst_id].append((start, end, ense_id))
                transsize[enst_id] += end - start

    for gene_id, thelist in list(manelist.items()):

        mane = (
            sorted(thelist[0], key=lambda x: transsize[x], reverse=True)
            + sorted(thelist[1], key=lambda x: transsize[x], reverse=True)
        )

        if not mane:
            warning(f"gene ID '{gene_id}' has no transcript/exon record in annotation")
            continue

        manelist[gene_id] = [mane[0], transtoexons[mane[0]]]

    return genetoanno, transtoexons, transtogenes, manelist, exon_ids


def readtable(inf):

    records = cl.defaultdict(list)
    with open_text(inf) as fin:
        for ln, line in enumerate(fin, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                warning(f"line {ln} has fewer than 5 columns; skipping")
                continue

            c1 = parts[0]
            c4 = [x for x in parts[4].split(";") if len(x)]

            # Split 5th column on ';' and parse each non-empty token
            for token in c4:
                try:
                    field1, field2, field3 = token.split(":")
                except ValueError:
                    warning(f"line {ln} token '{token}' is not formatted as ID:field:value; skipping")
                    continue

                note = NOTE_PATTERN.search(field1)
                note = f"({note.group(1)})" if note else ""
                field1 = field1.split("(")[0]
                records[field1].append((note, field2, float(field3), c1))

    return records


def summary(records, genetoanno, transtoexons, transtogenes, manelist, exon_ids, expected_gencode_version):

    gene_records = cl.defaultdict(list)
    gene_cn = cl.defaultdict(float)

    known_gene_ids = set(genetoanno)
    known_transcript_ids = set(transtogenes)
    known_exon_ids = set(exon_ids)
    warned_ids = set()
    check_msg = gencode_check_message(expected_gencode_version)

    for trans, data in records.items():

        if trans not in transtogenes:
            if trans not in warned_ids:
                if trans in known_gene_ids:
                    warning(
                        f"gene ID '{trans}' was found in the annotation record, "
                        f"but this script summarizes transcript IDs; skipping; {check_msg}"
                    )
                elif trans in known_exon_ids:
                    warning(
                        f"exon ID '{trans}' was found in the annotation record, "
                        f"but this script summarizes transcript IDs; skipping; {check_msg}"
                    )
                elif trans in known_transcript_ids:
                    warning(
                        f"transcript ID '{trans}' was found in the annotation record, "
                        f"but no gene/exon mapping was available; skipping; {check_msg}"
                    )
                else:
                    warning(
                        f"gene/transcript/exon ID '{trans}' was not found in annotation record; "
                        f"skipping; {check_msg}"
                    )
                warned_ids.add(trans)
            continue

        genename = transtogenes[trans]
        exons = transtoexons[trans]

        for (note, field2, field3, c1) in data:

            if float(field3) < 0.95:
                continue

            if note == "":
                mane = trans
                mane_exons = exons
                exon_start = 1
                exon_end = len(exons)
            else:
                if genename not in manelist:
                    warning(
                        f"gene ID '{genename}' for transcript ID '{trans}' was not found in "
                        f"MANE/transcript record; skipping; {check_msg}"
                    )
                    continue
                mane, mane_exons = manelist[genename]
                exon_start, exon_end = [int(x) for x in note[1:-1].split("-")]

            genesize = sum([x[1] - x[0] for x in mane_exons])

            if exon_end > len(exons):
                warning(
                    f"exon range '{exon_start}-{exon_end}' for transcript ID '{trans}' "
                    f"exceeds available exon records ({len(exons)}); skipping; {check_msg}"
                )
                continue

            exon_range = range(exon_start - 1, exon_end)
            exon_coordi = [exons[i] for i in exon_range]
            exon_overlaps, exon_index = interval_intersection(mane_exons, exons)

            if len(exon_index) == 0:
                warning(f"no exon overlap found for transcript ID '{trans}' in annotation record; skipping; {check_msg}")
                continue

            exon_index = "_".join(compress_ranges(sorted(list(exon_index))))
            totalsize = sum([x[1] - x[0] for x in exon_overlaps])
            cn = totalsize / max(1, genesize)
            gene_cn[mane] += cn

            gene_records[mane].append((c1, genesize, len(mane_exons), f"Exon_Index:{exon_index}", field3))

    return gene_cn, gene_records


def run(args):

    metadata = detect_gencode_metadata(args.gene)
    validate_gencode_metadata(args, metadata)

    records = readtable(args.input)

    genetoanno, transtoexons, transtogenes, manelist, exon_ids = readgff3(args.gene, records)

    gene_cn, gene_records = summary(
        records,
        genetoanno,
        transtoexons,
        transtogenes,
        manelist,
        exon_ids,
        str(args.expected_gencode_version),
    )

    # Write results
    results = []

    for trans, cn in gene_cn.items():

        genename = transtogenes[trans]
        genesize, exonnum = gene_records[trans][0][1:3]

        record = sorted(
            gene_records[trans],
            key=lambda x: (
                int(x[3].split(":")[1].split("-")[0]),
                int(x[3].split(":")[1].split("-")[-1]),
                -x[2],
                x[0],
            ),
        )

        record_str = ";".join([f"{x[0]}:{x[3]}:{x[4]}" for x in record])
        results.append(
            f"{trans}\t{genename}\t{genetoanno[genename][0]}\t{genetoanno[genename][1]}\t"
            f"{genesize}\t{exonnum}\t{cn}\t{record_str}\n"
        )

    results = sorted(results, key=lambda x: (x.split("\t")[2], x))

    with open(args.output, "w", encoding="utf-8") as fout:
        fout.write("genename\tid\ttype\ttotal_size\ttotal_exons\taggregate_copy_number\talleles\n")

        for line in results:
            fout.write(line)


def main():
    parser = argparse.ArgumentParser(
        description="Getting aggregate copy number calling."
    )
    parser.add_argument("-i", "--input", required=True, help="Input text file")
    parser.add_argument("-g", "--gene", required=True, help="Input genecode database")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    parser.add_argument(
        "--expected-genome-build",
        default="hg38",
        help="Expected genome build for the GENCODE annotation. Default: hg38",
    )
    parser.add_argument(
        "--expected-gencode-version",
        default="37",
        help="Expected/default GENCODE version. Default: 37",
    )
    parser.add_argument(
        "--strict-gencode-version",
        action="store_true",
        help="Stop if the detected GENCODE version or genome build does not match the expected values.",
    )
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
