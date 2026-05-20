import pysam
import sys
import csv

def process_bam(input_bam):
    # Open the BAM file for reading
    samfile = pysam.AlignmentFile(input_bam, "rb")
    output_file = input_bam.replace(".bam", ".paired.csv")

    # Define the CSV columns
    header = [
        "QNAME", "R1_FLAG", "R1_POS", "R1_CIGAR", "R1_MAPQ", "R1_NM", "R1_AS",
        "R2_FLAG", "R2_POS", "R2_CIGAR", "R2_MAPQ", "R2_NM", "R2_AS", "INSERT_SIZE"
    ]

    # Dictionary to hold the first read of a pair
    cache = {}

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)

        for read in samfile:
            # Skip unmapped or non-paired reads if necessary
            if not read.is_paired:
                continue

            qname = read.query_name

            if qname in cache:
                # We found the partner!
                r1 = cache.pop(qname)
                r2 = read

                # Ensure we know which is Read 1 and which is Read 2
                # (If the first one we saw was actually R2, we swap them)
                if r2.is_read1:
                    r1, r2 = r2, r1

                # Extract optional tags safely using .get_tag()
                # We use .get() style defaults to avoid errors if tags are missing
                def get_val(r, tag):
                    try:
                        return r.get_tag(tag)
                    except KeyError:
                        return "-"

                row = [
                    qname,
                    r1.flag, r1.reference_start + 1, r1.cigarstring, r1.mapping_quality, get_val(r1, "NM"), get_val(r1, "AS"),
                    r2.flag, r2.reference_start + 1, r2.cigarstring, r2.mapping_quality, get_val(r2, "NM"), get_val(r2, "AS"),
                    abs(r1.template_length)
                ]
                writer.writerow(row)
            else:
                # Store the current read and wait for its mate
                cache[qname] = read

    samfile.close()
    print(f"Done! Created {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python bam_to_paired_csv.py <input.bam>")
    else:
        process_bam(sys.argv[1])
