#!/usr/bin/env python3
"""Mask positions in every sequence of a FASTA alignment using a BED file.

The chromosome/name column of the BED file is ignored, so the same intervals
are applied to all sequences regardless of their names. Alignment column i is
assumed to correspond to 0-based genome position i (BED coordinates are
0-based, half-open).

Example:
    python mask_alignment.py alignment.fasta mask.bed alignment.masked.fasta
"""

import argparse
import sys


def read_bed(path):
    """Return a list of (start, end) tuples from a BED file (name column ignored)."""
    intervals = []
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"{path}:{lineno}: expected at least 3 columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                raise ValueError(f"{path}:{lineno}: start/end must be integers")
            if start < 0 or end < start:
                raise ValueError(f"{path}:{lineno}: invalid interval {start}-{end}")
            intervals.append((start, end))
    return intervals


def read_fasta(path):
    """Yield (header, sequence) for each record in a FASTA file."""
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:], []
            elif line.strip():
                if name is None:
                    raise ValueError(f"{path}: sequence data found before first header")
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks)


def write_fasta_record(fout, name, seq, width):
    fout.write(f">{name}\n")
    if width and width > 0:
        for i in range(0, len(seq), width):
            fout.write(seq[i:i + width] + "\n")
    else:
        fout.write(seq + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Mask positions from a BED file in all sequences of a FASTA alignment."
    )
    parser.add_argument("input_fasta", help="input alignment FASTA")
    parser.add_argument("bed", help="BED file of regions to mask (chrom column ignored)")
    parser.add_argument("output_fasta", help="output masked FASTA")
    parser.add_argument("-c", "--mask-char", default="N",
                        help="character used for masking (default: N)")
    parser.add_argument("-w", "--width", type=int, default=0,
                        help="line width for output sequences (default: 0 = single line)")
    args = parser.parse_args()

    if len(args.mask_char) != 1:
        parser.error("--mask-char must be a single character")

    intervals = read_bed(args.bed)
    max_end = max((end for _, end in intervals), default=0)

    aln_len = None
    n_records = 0

    with open(args.output_fasta, "w") as fout:
        for name, seq in read_fasta(args.input_fasta):
            if aln_len is None:
                aln_len = len(seq)
                if max_end > aln_len:
                    sys.exit(f"Error: BED extends to {max_end}, but alignment length is {aln_len}")
            elif len(seq) != aln_len:
                sys.exit(f"Error: sequence '{name}' has length {len(seq)}, "
                         f"expected {aln_len} (is this an alignment?)")

            masked = bytearray(seq.encode())
            fill = args.mask_char.encode()
            for start, end in intervals:
                masked[start:end] = fill * (end - start)

            write_fasta_record(fout, name, masked.decode(), args.width)
            n_records += 1

    if n_records == 0:
        sys.exit("Error: no sequences found in input FASTA")

    print(f"Masked {len(intervals)} interval(s) in {n_records} sequence(s) "
          f"(alignment length {aln_len}) -> {args.output_fasta}", file=sys.stderr)


if __name__ == "__main__":
    main()
