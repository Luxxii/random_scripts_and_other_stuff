#!/bin/env python

import argparse
import csv
import sys
from ast import literal_eval

import pandas as pd
import protgraph as pg
import tqdm

csv.field_size_limit(sys.maxsize)

HYDROGEN_MONO_MASS = 1.007825035
NEUTRON_MONO_MASS = 1.00866491588
OXYGEN_MONO_MASS = 15.994915
AMINO_ACID_DICT = pg.aa_masses_annotation._get_mass_dict(factor=1, type=float)

def parse_args():
    """ Argument parsing """
    parser = argparse.ArgumentParser()
    parser.add_argument("-input_sequences", help="Sequences (like SpikeIns) in csv, which should be matched to the found features in UnbeQuant. This should contain the following column_headers: Accession, Gene, Sequence")
    parser.add_argument("-input_unbequant", help="UnbeQuant Input (in tsv), which contains the column mz_start and mz_end for the features/isotopes")
    parser.add_argument("-fixmod", help="Fixed Modifications -> C:57.XXX.  Apply this multiple times to add more modifications", action="append", default=[])
    parser.add_argument("-varmod", help="Variable Modifications -> M:18.XXX.  Apply .this multiple times to add more modifications", action="append", default=[])
    parser.add_argument("-max_charge", help="Maximum Charge to consider")
    parser.add_argument("-chech_x_first_isotopes", help="match only to the first X isotopes", default=3)
    parser.add_argument("-out_tsv", help="Output-tsv-file, summarizing features")
    parser.add_argument("-only_output_hits", help="Set if only hits should be saved in the final output", action="store_true", default=False)
    return parser.parse_args()

def get_masses(in_file, varmod, fixmod, max_charge, max_first_isotopes):
    """ Get the masses of the sequences in the input file, including variable/fixed modifications and up to n isotopes. """
    with open(in_file, "r") as in_filtered:
        csv_in = csv.reader(in_filtered, delimiter="\t")

        filtered_header = next(csv_in)
        accession_idx = filtered_header.index("Accession")
        gene_idx = filtered_header.index("Gene")
        sequence_idx = filtered_header.index("Sequence")

        entries = dict()
        for l in csv_in:
            # Calcualte Mass and FIXMOD
            mono_mass = sum([AMINO_ACID_DICT[x][0] for x in l[sequence_idx]]) + 2*HYDROGEN_MONO_MASS + 1*OXYGEN_MONO_MASS # Mass + H2O
            fm_mass_shift = 0
            for fm in fixmod:
                fix_count = l[sequence_idx].count(fm.split(":")[0])
                fm_mass_shift += fix_count*float(fm.split(":")[1])
            mono_mass += fm_mass_shift  # Include FIXMOD
        
            # Get all possible VARMOD and their combinations
            varmod_possibilities = [[     ((i)*float(vm.split(":")[1]), str(i) + "x" + vm if i != 0 else "") for i in range(l[sequence_idx].count(vm.split(":")[0])+1)    ]  for vm in varmod]
            
            # Now combine every two lists, with each other to get all combinations of VARMODs
            while len(varmod_possibilities) > 1:
                new_varmod_possibilities = []
                for x in varmod_possibilities[0]:
                    for y in varmod_possibilities[1]:
                        new_varmod_possibilities.append((x[0] + y[0], x[1] + ";" + y[1]))
                varmod_possibilities = [new_varmod_possibilities] + varmod_possibilities[2:]
            
            # Then iterate through it to get all distinct mono masses
            mono_masses_varmod = []
            mono_masses_applied_mods = []
            for delta, mod in varmod_possibilities[0]:
                mono_masses_varmod.append(mono_mass + delta)
                mono_masses_applied_mods.append(mod)
            # Included VARMODs

            # Get all the single charges:
            mono_masses_charged = []
            mono_masses_charge_state = []
            mono_masses_applied_mods_and_charge_state = []
            for c in range(1, int(max_charge) + 1):
                for mmv, mmvam in zip(mono_masses_varmod, mono_masses_applied_mods):
                    for i in range(0, int(max_first_isotopes)):  # We check for the first x isotopes
                        mmv_with_iso = mmv + i*(NEUTRON_MONO_MASS)
                        mono_masses_charged.append(
                            (c*HYDROGEN_MONO_MASS + mmv_with_iso) / c
                        )
                        mono_masses_charge_state.append(c)
                        mono_masses_applied_mods_and_charge_state.append(mmvam)

            # Add the entries
            entries[l[sequence_idx]] = (
                l[gene_idx],  # Gene Information
                l[accession_idx], # Accession Information
                mono_masses_charged, # List of possible weights
                mono_masses_charge_state, # List of possible charges
                mono_masses_applied_mods_and_charge_state # List of applied modifications
            )

        return entries


if __name__ == "__main__":
    args = parse_args()

    # First load possible masses from table
    seq_mass_dict = get_masses(
        args.input_sequences,
        args.varmod,
        args.fixmod,
        args.max_charge,
        args.chech_x_first_isotopes
    )

    # Create a pandas DataFrame for quicker accession to all masses
    pd_dict_frame = {"gene": [], "accession": [], "sequence": [], "masses": [], "charges": [], "modifications": []}
    for seq, val in seq_mass_dict.items():
        pd_dict_frame["gene"].extend([val[0]]*len(val[2]))
        pd_dict_frame["accession"].extend([val[1]]*len(val[2]))
        pd_dict_frame["sequence"].extend([seq]*len(val[2]))
        pd_dict_frame["masses"].extend(val[2])
        pd_dict_frame["charges"].extend(val[3])
        pd_dict_frame["modifications"].extend(val[4])
    df_seq_mass = pd.DataFrame(pd_dict_frame)

    with open(args.input_unbequant, "r") as in_anon, open(args.out_tsv, "w") as out_file:
        csv_in = csv.reader(in_anon, delimiter="\t")
        csv_out = csv.writer(out_file, delimiter="\t")

        headers = next(csv_in)
        csv_out.writerow(headers + ["hitsearch_accession", "hitsearch_gene", "hitsearch_sequence", "hitesearch_modifications"])

        # Get the indices of the needed columns
        id_idx = headers.index("openms_ceid")
        mz_start_idcs = [x for x in headers if x.endswith("_____l_mz_start")]
        mz_start_idcs = [headers.index(x) for x in mz_start_idcs]
        mz_end_idcs = [x for x in headers if x.endswith("_____l_mz_end")]
        mz_end_idcs = [headers.index(x) for x in mz_end_idcs]
        charge_idcs = [x for x in headers if x.endswith("_____charge")]
        charge_idcs = [headers.index(x) for x in charge_idcs]

        # Get for each isotope the max and minimum interval
        all_features = dict()
        for l in tqdm.tqdm(csv_in, unit="lines"):
            starts = [literal_eval(l[mzs]) for mzs in mz_start_idcs if l[mzs]]
            ends = [literal_eval(l[mze]) for mze in mz_end_idcs if l[mze]]

            start_end = [[float("infinity"), -float("infinity")]]
            for s, e in zip(starts, ends):
                while len(s) > len(start_end):
                    start_end.append([float("infinity"), -float("infinity")])
                for idx, (iso_s, iso_e) in enumerate(zip(s,e)):
                    iso_s, iso_e = float(iso_s), float(iso_e)
                    if iso_s < start_end[idx][0]:
                        start_end[idx][0] = iso_s
                    if iso_e > start_end[idx][1]:
                        start_end[idx][1] = iso_e

            # Get the matching of a peptide to be searched to a feature
            hitsearch_accession = []
            hitsearch_gene = []
            hitsearch_sequence = []
            hitesearch_modifications = []
            min_mz, max_mz = start_end[0] # We only look at the lightest isotope
            mz_charges = list(set([l[x] for x in charge_idcs]))

            # Get all matching peptides from spikeins
            for mz_c in mz_charges:
                if mz_c != "":
                    c_frame = df_seq_mass[df_seq_mass["charges"] == int(mz_c)]
                    c_mz_frame = c_frame[c_frame["masses"].between(min_mz, max_mz)]

                    if not c_mz_frame.empty:
                        hitsearch_accession.extend(c_mz_frame["accession"].tolist())
                        hitsearch_gene.extend(c_mz_frame["gene"].tolist())
                        hitsearch_sequence.extend(c_mz_frame["sequence"].tolist())
                        hitesearch_modifications.extend(c_mz_frame["modifications"].tolist())

            # Write all the matches of the current line into the final output:
            if args.only_output_hits:
                if len(hitsearch_accession) > 0:
                    csv_out.writerow(
                        l + [
                            "|".join(hitsearch_accession), "|".join(hitsearch_gene), "|".join(hitsearch_sequence), "|".join(hitesearch_modifications)
                        ]
                    )
            else:
                csv_out.writerow(
                    l + [
                        "|".join(hitsearch_accession), "|".join(hitsearch_gene), "|".join(hitsearch_sequence), "|".join(hitesearch_modifications)
                    ]
                )
