import argparse
import csv
import sys

import tqdm
import numpy as np
from scipy.stats import ttest_ind
from scipy.stats import false_discovery_control
import pandas as pd
from ast import literal_eval

csv.field_size_limit(sys.maxsize)

NONE_VALUES = {
    "nan",
    "inf",
    "#n/a",
    "n/a",
    "na",
    "",
}


def parse_ms2_scans(cells: str, origin_str: str):
    """ Parse the MS2 scans from a list of cell entries """
    ms2s = set()
    ms2_with_origin = []
    for cell, origin in zip(cells, origin_str):
        if cell:
            if cell.lower() in NONE_VALUES:
                ms2_with_origin.append(None)
                continue
            else:
                try:
                    ms2s = ",".join(str(x) for x in literal_eval(cell))
                    if ms2s:
                        ms2_with_origin.append(origin + ":" + ms2s)
                except (ValueError, SyntaxError):
                    # if the cell is not parsable it is a malformed list string element (we skip it)
                    ms2_with_origin.append(None)
    return ms2_with_origin


def parse_pep_prot_ident(cells: str):
    """ Parse the peptide/protein identification from a list of cell entries """
    peps_prots = set()
    for cell in cells:
        if cell:
            if cell.lower() in NONE_VALUES:
                continue
            else:
                try:
                    for c in literal_eval(cell):
                        peps_prots.add(c)
                except (ValueError, SyntaxError):
                    # if the cell is not parsable it is a malformed list string element (we skip it)
                    continue

    return peps_prots

def parse_charge(cell: str):
    """ Parse the charge value from a cell, returns None if not parsable"""
    if cell:
        try:
            return int(cell)
        except ValueError:
            # if the cell is not parsable (e.g. #N/A or similar, we return None)
            return None
    else:
        return None

def parse_intensity(cell: str):
    """ Parse the intensity value from a cell, returns None if the cell is empty or not parsable """
    if cell:
        if cell.lower() == "nan":
            return None
        if cell.lower() == "inf":
            return None
        try:
            return float(cell)
        except ValueError:
            # if the cell is not parsable (e.g. #N/A or similar, we return None)
            return None
    else:
        return None


def parse_args():
    """ Parse Arguments """
    parser = argparse.ArgumentParser(
        description="A small script to calculate the p-value in a CSV-file across two groups (same size), using column-headers."
    )

    # Input CSV
    parser.add_argument(
        "--input_unbequant_tsv", "-in", type=argparse.FileType("r"),
        help="The Output of Unbequant (intensities can be normalized before)."
    )

    parser.add_argument(
        "--groupA", "-gA", type=str, action="append",
        help="Group A. Specify here the column-headers, which belong to this group"
    )
    parser.add_argument(
        "--groupB", "-gB", type=str, action="append",
        help="Group B. Specify here the column-headers, which belong to this group"
    )

    parser.add_argument(
        "--percentage_of_missingness_per_group_allowed", "-p", type=float, default=0,
        help="A percentage (0 - 1), which is allowed to be missing in a group. Set it to 0.3 to allow 30% of missing values in a group"
    )

    parser.add_argument(
        "--log2_intensities_provided", "-lip", type=bool, default=True,
        help="Are there log2 intensities provided? If yes, set this to 1. "
    )

    # Output-CSV
    parser.add_argument(
        "--output", "-out", type=str, default=sys.stdout,
        help="Output-TSV containing the same columns and rows, with additional p_value, corrected p_value and fold_change"
        " Defaults to stdout."
    )

    # Output-CSV
    parser.add_argument(
        "--output_volcano", "-volcano", type=str, default=None,
        help="Output-HTML of the volcano plot. "
        " Defaults to None (--> No plots generated)."
    )

    return parser.parse_args()



def coloring_volcano(x):
    if x == 6:
        # is not identified and not significant
        return "Not Identified and not Significant"
    elif x == 7:
        # is identified and not significant
        return "Identified and not Significant"
    elif x == 11:
        # is not identified and significant
        return  "Not Identified and Significant"
    elif x == 12:
        # is identified and significant
        return "Identified and Significant"



if __name__ == "__main__":
    # Parse arguments
    args = parse_args()

    # TODO should they be in the same size ?
    # if args.groupA is None or args.groupB is None or \
    #     len(args.groupA) != len(args.groupB) or len(args.groupA) == 0:
    #     raise Exception("Groups need to be specified and need to be of same length")

    # Initialize csv reader/writer
    tsv_in = csv.reader(args.input_unbequant_tsv, delimiter="\t")
    # Get indices of the columns
    header = next(tsv_in)
    gA_idx = [header.index(x) for x in args.groupA]
    gB_idx = [header.index(x) for x in args.groupB]
    openmsid_idx = header.index("openms_ceid")
    pep_ident_idx = [header.index(x) for x in header if x.endswith("_____l_pep_ident")]
    raw_pep_ident_idx = [header.index(x) for x in header if x.endswith("_____l_raw_pep_ident")]
    prot_ident_idx = [header.index(x) for x in header if x.endswith("_____l_prot_ident")]
    charge_idx = [header.index(x) for x in header if x.endswith("_____charge")]
    ms2_scans_idx = [header.index(x) for x in header if x.endswith("_____l_ms2_scans")]
    ms2_scans_origin = [x[:-len("_____l_ms2_scans")] for x in header if x.endswith("_____l_ms2_scans")]
    first_iso_global_min_mz_idx = header.index("first_iso_global_min_mz")
    first_iso_global_max_mz_idx = header.index("first_iso_global_max_mz")
    first_iso_global_min_rt_idx = header.index("first_iso_global_min_rt")
    first_iso_global_max_rt_idx = header.index("first_iso_global_max_rt")
    feature_global_min_mz_idx = header.index("feature_global_min_mz")
    feature_global_max_mz_idx = header.index("feature_global_max_mz")
    feature_global_min_rt_idx = header.index("feature_global_min_rt")
    feature_global_max_rt_idx = header.index("feature_global_max_rt")
    # Get additional headers, which were addded after the UnbeQuant analysis
    non_unbequant_headers = set([header.index(x) for x in header if "_____" not in x]) - set([openmsid_idx, first_iso_global_min_mz_idx, first_iso_global_max_mz_idx, first_iso_global_min_rt_idx, first_iso_global_max_rt_idx, feature_global_min_mz_idx, feature_global_max_mz_idx, feature_global_min_rt_idx, feature_global_max_rt_idx])

    # Iterate over each line and generate pd Frame
    dict_for_pd = {**dict(
        openmsid=[],
        ttest_ind_pvalue=[],
        ttest_ind_statistic=[],
        fold_change_A_div_B=[],
        missing_values_in_a=[],
        missing_values_in_b=[],
        charge=[],
        first_iso_global_min_mz=[],
        first_iso_global_max_mz=[],
        first_iso_global_min_rt=[],
        first_iso_global_max_rt=[],
        feature_global_min_mz=[],
        feature_global_max_mz=[],
        feature_global_min_rt=[],
        feature_global_max_rt=[],
        pep_ident=[],
        raw_pep_ident=[],
        prot_ident=[],
        ms2s=[],
    ), **{header[x]: [] for x in non_unbequant_headers}}
    additional_headers = [header[x] for x in non_unbequant_headers]
    for line in tqdm.tqdm(tsv_in):

        # Get Cell Contents
        ga = [parse_intensity(line[x]) for x in gA_idx]
        gb = [parse_intensity(line[x]) for x in gB_idx]

        # Check if we have enough values
        if ga.count(None)/len(ga) > args.percentage_of_missingness_per_group_allowed:
            # Too many missing values in group A
            continue
        if gb.count(None)/len(gb) > args.percentage_of_missingness_per_group_allowed:
            # Too many missing values in group B
            continue

        missing_ga = ga.count(None)/len(ga)
        missing_gb = gb.count(None)/len(gb)

        ga = [x for x in ga if x is not None]
        gb = [x for x in gb if x is not None]

        try:
            # Write pvalue and statistic
            res = ttest_ind(ga, gb)
            # Add to the dictionary missing information and other static stuff
            dict_for_pd["openmsid"].append(line[openmsid_idx])
            dict_for_pd["pep_ident"].append(";".join(parse_pep_prot_ident([line[x] for x in pep_ident_idx])))
            dict_for_pd["raw_pep_ident"].append(";".join(parse_pep_prot_ident([line[x] for x in raw_pep_ident_idx])))
            dict_for_pd["prot_ident"].append(";".join(parse_pep_prot_ident([line[x] for x in prot_ident_idx])))
            dict_for_pd["charge"].append((set([parse_charge(line[x]) for x in charge_idx]) - {None}).pop())
            dict_for_pd["ms2s"].append(";".join(parse_ms2_scans([line[x] for x in ms2_scans_idx],  ms2_scans_origin)))
            dict_for_pd["ttest_ind_pvalue"].append(res.pvalue)
            dict_for_pd["ttest_ind_statistic"].append(res.statistic)
            if args.log2_intensities_provided:
                dict_for_pd["fold_change_A_div_B"].append(np.log2(np.mean([2**x for x in ga]) / np.mean([2**x for x in gb])))
            else:
                dict_for_pd["fold_change_A_div_B"].append(np.log2(np.mean(ga) / np.mean(gb)))
            dict_for_pd["missing_values_in_a"].append(missing_ga)
            dict_for_pd["missing_values_in_b"].append(missing_gb)
            dict_for_pd["first_iso_global_min_mz"].append(float(line[first_iso_global_min_mz_idx]))
            dict_for_pd["first_iso_global_max_mz"].append(float(line[first_iso_global_max_mz_idx]))
            dict_for_pd["first_iso_global_min_rt"].append(float(line[first_iso_global_min_rt_idx]))
            dict_for_pd["first_iso_global_max_rt"].append(float(line[first_iso_global_max_rt_idx]))
            dict_for_pd["feature_global_min_mz"].append(float(line[feature_global_min_mz_idx]))
            dict_for_pd["feature_global_max_mz"].append(float(line[feature_global_max_mz_idx]))
            dict_for_pd["feature_global_min_rt"].append(float(line[feature_global_min_rt_idx]))
            dict_for_pd["feature_global_max_rt"].append(float(line[feature_global_max_rt_idx]))
            for h in non_unbequant_headers:
                dict_for_pd[header[h]].append(line[h])
        except Exception as e:
            # Could not calculate write nones
            print("Could not calculate ttest_ind for rowkey: ", line[openmsid_idx], "Error: ", e)


    # Generate pandas dataframe
    df = pd.DataFrame(dict_for_pd)

    # Sort by ascending p value:
    df = df.sort_values(by=["ttest_ind_pvalue"], ascending=True)

    # Do p-value correction:
    df["ttest_ind_corrected_pvalue"] = false_discovery_control(df["ttest_ind_pvalue"], method="bh")

    # Calculate the p-value for volcano
    df["-log10(pvalue)"] = -np.log10(df["ttest_ind_pvalue"])

    # Save the dataframe
    df.to_csv(args.output, sep="\t", index=False)

    # Save the volcano plot
    if args.output_volcano:
        df["pep_ident_hovers"] = df["pep_ident"].apply(lambda x: x[:30] + (x[30:] and ".."))
        df["prot_ident_hovers"] = df["prot_ident"].apply(lambda x: x[:30] + (x[30:] and ".."))
        df["raw_pep_ident_hovers"] = df["raw_pep_ident"].apply(lambda x: x[:30] + (x[30:] and ".."))
        df["ms2s_hovers"] = df["ms2s"].apply(lambda x: x[:30] + (x[30:] and ".."))
        for h in additional_headers:
            df[h + "_hovers"] = df[h].apply(lambda x: str(x)[:30] + (str(x)[30:] and ".."))
        df["is_identified"] = df["pep_ident"].apply(lambda x: 2 if x != "" else 1)
        df["is_significat"] = (df["ttest_ind_corrected_pvalue"] < 0.05) * 5 + 5
        df["is_ident_significat"] = df["is_identified"] + df["is_significat"]
        df["color"] = df["is_ident_significat"].apply(coloring_volcano)
        color_discrete_map = {
            'Not Identified and not Significant': '#365575',
            'Identified and not Significant"': '#4f7832',
            'Not Identified and Significant': '#076ddb',
            'Identified and Significant': '#59cc06'
        }

        import plotly.express as px
        fig = px.scatter(df, x="fold_change_A_div_B", y="-log10(pvalue)", color="color", color_discrete_map=color_discrete_map, hover_data={**{
            "color" : False,
            "openmsid" : True,
            "ttest_ind_pvalue" : True,
            "ttest_ind_statistic" : True,
            "fold_change_A_div_B" : True,
            "missing_values_in_a" : True,
            "missing_values_in_b" : True,
            "charge" : True,
            "first_iso_global_min_mz" : True,
            "first_iso_global_max_mz" : True,
            "first_iso_global_min_rt" : True,
            "first_iso_global_max_rt" : True,
            "feature_global_min_mz" : True,
            "feature_global_max_mz" : True,
            "feature_global_min_rt" : True,
            "feature_global_max_rt" : True,
            "pep_ident_hovers" : True,
            "raw_pep_ident_hovers" : True,
            "prot_ident_hovers" : True,
            "ms2s_hovers" : True,
        } , **{h + "_hovers": True for h in additional_headers}})
        fig.add_hline(y=-np.log10(0.05), line_width=1, line_dash='dash')
        fig.add_vline(x=np.log2(2), line_width=1, line_dash='dash')
        fig.add_vline(x=np.log2(0.5), line_width=1, line_dash='dash')

        # Write volcano
        fig.write_html(args.output_volcano)
