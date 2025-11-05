#!/usr/bin/env python3
import argparse
import subprocess
import pandas as pd
import os
import tempfile


def ensure_repeatmasker_bed(rmsk_bed_path: str, rmsk_txt_gz: str = "rmsk.hg38.txt.gz"):
    """Check if RepeatMasker BED exists; if not, download and generate it."""
    if os.path.exists(rmsk_bed_path):
        print(f"RepeatMasker BED file found at {rmsk_bed_path}")
        return

    print("RepeatMasker BED not found. Downloading and preparing annotation...")

    rmsk_bed_tmp = rmsk_bed_path

    # Download RepeatMasker table
    subprocess.run(
        [
            "wget",
            "-O",
            rmsk_txt_gz,
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz",
        ],
        check=True,
    )

    # Convert to BED format using awk
    awk_cmd = (
        f"zcat {rmsk_txt_gz} | "
        "awk 'BEGIN{OFS=\"\\t\"} {print $6, $7, $8, $11, $2, $10}' "
        f"> {rmsk_bed_tmp}"
    )
    subprocess.run(awk_cmd, shell=True, check=True)

    # Remove compressed text file after conversion
    os.remove(rmsk_txt_gz)
    print(f"RepeatMasker BED created at {rmsk_bed_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute repeat-overlap feature for E–G pairs using RepeatMasker annotations (hg38)."
    )
    parser.add_argument(
        "--e2g",
        required=True,
        help="Path to E–G pairs TSV/TSV.GZ file (must contain ElementChr, ElementStart, ElementEnd, ElementName)",
    )
    parser.add_argument(
        "--rmsk",
        required=False,
        default="repeatmasker_hg38.bed",
        help="Path to RepeatMasker BED file (will be downloaded if missing)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output file path for table with repeat overlap feature (TSV or TSV.GZ)",
    )
    args = parser.parse_args()

    e2g_file = args.e2g
    rmsk_file = args.rmsk
    output_file = args.output

    # Step 1. Ensure RepeatMasker BED is available
    rmsk_txt_gz = "rmsk.hg38.txt.gz"
    ensure_repeatmasker_bed(rmsk_file, rmsk_txt_gz)

    # Step 2. Prepare temporary BED file for unique elements
    print("Extracting and deduplicating genomic coordinates from E–G file...")
    tmp_bed = tempfile.NamedTemporaryFile(delete=False, suffix=".bed")
    tmp_bed_name = tmp_bed.name
    tmp_bed.close()

    e2g_df = pd.read_csv(e2g_file, sep="\t", compression="infer")
    bed_df = e2g_df[
        ["ElementChr", "ElementStart", "ElementEnd", "ElementName"]
    ].drop_duplicates(subset=["ElementChr", "ElementStart", "ElementEnd"])
    bed_df.to_csv(tmp_bed_name, sep="\t", header=False, index=False)

    # Step 3. Run bedtools coverage
    print("Running bedtools coverage...")
    tmp_cov = tempfile.NamedTemporaryFile(delete=False, suffix=".bed")
    tmp_cov_name = tmp_cov.name
    tmp_cov.close()

    bedtools_cmd = ["bedtools", "coverage", "-a", tmp_bed_name, "-b", rmsk_file]
    with open(tmp_cov_name, "w") as out_cov:
        subprocess.run(bedtools_cmd, stdout=out_cov, check=True)

    # Step 4. Load and process repeat overlap results
    print("Reading bedtools output and merging with E–G pairs...")
    cols = [
        "chr",
        "start",
        "end",
        "name",
        "length",
        "bases_covered",
        "fraction_covered",
    ]
    repeat_cov = pd.read_csv(tmp_cov_name, sep="\t", names=cols)

    repeat_cov["repeat_overlap_fraction"] = repeat_cov["fraction_covered"].fillna(0)
    repeat_cov["repeat_overlap_bp"] = repeat_cov["bases_covered"].fillna(0)
    repeat_cov["ElementName"] = (
        repeat_cov["chr"].astype(str)
        + ":"
        + repeat_cov["start"].astype(str)
        + "-"
        + repeat_cov["end"].astype(str)
    )

    # Step 5. Merge with original E–G table
    features = e2g_df.merge(
        repeat_cov[["ElementName", "repeat_overlap_fraction"]],
        on="ElementName",
        how="left",
    )
    features["repeat_overlap_fraction"].fillna(0, inplace=True)

    # Step 6. Write final output
    print(f"Writing output to {output_file}")
    features.to_csv(output_file, sep="\t", index=False, compression="infer")

    # Step 7. Cleanup
    os.remove(tmp_bed_name)
    os.remove(tmp_cov_name)
    print("Done. Temporary files removed.")


if __name__ == "__main__":
    main()
