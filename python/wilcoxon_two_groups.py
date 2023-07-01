import argparse
import csv
import sys

from scipy.stats import wilcoxon

csv.field_size_limit(sys.maxsize)


def parse_args():
    """ Parse Arguments """
    parser = argparse.ArgumentParser(
        description="A small script to calculate the p-value in a CSV-file across two groups (same size), using column-headers."
    )

    # Input CSV
    parser.add_argument(
        "--input_csv", "-in", type=argparse.FileType("r"),
        help="The CSV-file, where two groups are present."
    )

    parser.add_argument(
        "--groupA", "-gA", type=str, action="append",
        help="Group A. Specify here the column-headers, which belong to this group"
    )
    parser.add_argument(
        "--groupB", "-gB", type=str, action="append",
        help="Group B. Specify here the column-headers, which belong to this group"
    )

    # Output-CSV
    parser.add_argument(
        "--output", "-out", type=argparse.FileType("w"), default=sys.stdout,
        help="Output-CSV containing the same columns and rows, with additional p_value and statistics column added to the end. "
        " Defaults to stdout."
    )

    return parser.parse_args()


if __name__ == "__main__":
    # Parse arguments
    args = parse_args()

    if args.groupA is None or args.groupB is None or \
        len(args.groupA) != len(args.groupB) or len(args.groupA) == 0:
        raise Exception("Groups need to be specified and need to be of same length")

    # Initialize csv reader/writer
    csv_in = csv.reader(args.input_csv)
    csv_out = csv.writer(args.output)

    # Get indices of the columns
    header = next(csv_in)
    gA_idx = [header.index(x) for x in args.groupA]
    gB_idx = [header.index(x) for x in args.groupB]

    csv_out.writerow(header + ["wilcoxon_p_value", "wilcoxon_statistc"])

    # Iterate over each line
    for l_idx, line in enumerate(csv_in):
        
        ga = [float(line[x]) for x in gA_idx]
        gb = [float(line[x]) for x in gB_idx]
        try:
            # Write pvalue and statistic
            res = wilcoxon(ga, gb)
            csv_out.writerow(line + [res[1], res[0]])
        except:
            # Could not calculate write nones
            csv_out.writerow(line + [None, None])
