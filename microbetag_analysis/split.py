#!/usr/bin/env python3
import os
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path


AMINOACIDS = [
    "L-Tyrosine",
    "L-Phenylalanine",
    "L-Tryptophan",
    "L-Aspartic Acid",
    "L-Arginine",
    "L-Leucine",
    "L-Glutamic acid",
    "L-Glutamine",
    "L-Histidine",
    "L-Lysine",
    "L-Methionine",
    "L-Proline",
    "L-Threonine",
    "L-Valine",
    "L-Alanine",
    "L-Asparagine",
    "L-Isoleucine",
    "L-Serine"
]

VITAMINS = [
    "Thiamine"   ,  # B1 →
    "Riboflavin" ,  # B2 →
    "Niacinamide",  # B3 →
    "Nicotinic acid",
    "NAD",
    "Pantothenic acid",  # B5 →
    "Pyridoxamine"    ,  # B6 →
    "Pyridoxine",
    "a-Tocopheryl acetate"  # E →
]

OTHER = [
    "Adenine",
    "Guanosine",
    "Lactic acid",
    "Sulfoacetic acid",
    "Putrescine"
]

class SplitDataset:
    def __init__(
        self,
        abundance_file,
        genome_info,
        metadata_file,
        metabolites_file,
        categories,
        threshold   = 0.2,
        norm_metabo = None,
        outdir      = None,
        skip        = []
    ):

        self.threshold  = threshold
        self.categories = categories

        outdir     = outdir if outdir is not None else os.getcwd()
        if not os.path.isabs(outdir):
            script_dir  = os.path.dirname(os.path.abspath(__file__))
            self.outdir = os.path.join(script_dir, outdir)

        # ------------------------------------------------------------
        # Load abundance + taxonomy and build classification string
        # ------------------------------------------------------------
        abundance = pd.read_csv(abundance_file)

        tax = (
            pd.read_csv(genome_info)
            .loc[:, ["genome", "domain", "phylum", "class", "order", "family", "genus", "species"]]
        )

        tax["classification"] = tax.iloc[:, 1:].astype(str).agg(";".join, axis=1)

        abundance = abundance.merge(tax, on="genome", how="left")

        # ------------------------------------------------------------
        # Load metadata (weird format where metadata is transposed)
        # ------------------------------------------------------------
        metadata = pd.read_csv(metadata_file, header=None).T

        # First row is header
        metadata.columns = metadata.iloc[0]
        metadata = metadata.drop(index=0)

        # First metadata column is the new index
        metadata = metadata.set_index(metadata.columns[0])

        # Store sample names for later
        sample_names = metadata.columns.tolist()

        # ------------------------------------------------------------
        # Warn about missing samples in abundance
        # ------------------------------------------------------------
        missing = set(metadata.columns) - set(abundance.columns)
        if missing:
            warnings.warn(f"Missing in abundance: {missing}")

        # ------------------------------------------------------------
        # Load metabolite sheets (Digesta, Annotations)
        # ------------------------------------------------------------
        digesta     = pd.read_excel(metabolites_file, sheet_name="Digesta")
        annotations = pd.read_excel(metabolites_file, sheet_name="Annotations")

        # Map Feature_ID → Curated.ID
        mapping = annotations.set_index("Feature_ID")["Curated.ID"]
        digesta["Feature_ID"] = digesta["Feature_ID"].map(mapping)

        # Keep amino acids + vitamins
        allowed_ids = set(AMINOACIDS + VITAMINS + OTHER)
        digesta_filtered = digesta[digesta["Feature_ID"].isin(allowed_ids)]

        # Convert to long format (animal → features)
        dig_long = (
            digesta_filtered
            .set_index("Feature_ID")
            .T
            .rename_axis("animal")
            .reset_index()
        )
        if norm_metabo == "log1p":
            num_cols = dig_long.select_dtypes(include=[np.number]).columns
            dig_long_norm = dig_long.copy()
            dig_long_norm[num_cols] = np.log1p(dig_long[num_cols])

        # ------------------------------------------------------------
        # Attach metabolite info to metadata
        # ------------------------------------------------------------
        metadata_T = metadata.T.copy()

        metadata_T["sample"]   = metadata_T.index
        metadata_T["animal_m"] = metadata_T["animal"] + "M"

        metadata_T = metadata_T.merge(dig_long, left_on="animal_m", right_on="animal", how="left")

        # ------------------------------------------------------------
        # Reorder abundance columns to match metadata
        # ------------------------------------------------------------
        ordered_samples = abundance.columns.intersection(metadata.columns)
        abundance       = abundance[["genome"] + ordered_samples.tolist() + ["classification"]]

        # ------------------------------------------------------------
        # Final transpose of metadata + restore original sample names
        # ------------------------------------------------------------
        metadata_final         = metadata_T.T
        metadata_final.columns = sample_names

        metadata_f = metadata_final[ordered_samples]
        # metadata_f = metadata_f.drop(index=["batch", "tissue", "animal_m", "animal_y", "sample"], errors="ignore")
        skip_rows = ["batch", "tissue", "animal_m", "animal_y", "sample"] + skip
        metadata_f = metadata_f.drop(index=skip_rows, errors="ignore")

        # Assign outputs
        self.metadata  = metadata_f
        self.abundance = abundance

    def split(self):

        for type in self.categories:

            # Get the treatment row as a Series
            cases = self.metadata.loc[type, :]

            for case in cases.unique():

                # Get column names (samples) matching the current treatment
                selected_samples = cases[cases == case].index

                # Subset the metadata to only those samples (i.e., columns)
                tr_metadata = self.metadata[selected_samples]

                # Get same samples from the abundance table
                tr_abd = self.abundance[["genome"] + tr_metadata.columns.tolist() + ["classification"]]

                df_filtered = SplitDataset.apply_prevalence(tr_abd, self.threshold)

                # Export those subset of the dataframes to files
                abd_outfile  = "_".join(["abd_prev", str(self.threshold), str(type), str(case)])
                abd_outfile += ".csv"
                meta_outfile = "_".join(["metadata_prev", str(self.threshold), str(type), str(case)])
                meta_outfile += ".tsv"

                if self.outdir:
                    os.makedirs(self.outdir, exist_ok=True)
                df_filtered.to_csv(Path(self.outdir) / abd_outfile, index=False)
                tr_metadata.to_csv(Path(self.outdir) / meta_outfile, index_label="sample", sep="\t")

    def overall(self):

        # Export those subset of the dataframes to files
        abd_outfile  = "_".join(["abd_prev", str(self.threshold), "overall"])
        abd_outfile += ".csv"
        meta_outfile = "_".join(["metadata_prev", str(self.threshold), "overall"])
        meta_outfile += ".tsv"

        df_filtered = SplitDataset.apply_prevalence(self.abundance, self.threshold)

        if self.outdir:
            os.makedirs(self.outdir, exist_ok=True)
        df_filtered.to_csv(Path(self.outdir) / abd_outfile, index=False)
        self.metadata.to_csv(Path(self.outdir) / meta_outfile, index_label="sample", sep="\t")

    @staticmethod
    def apply_prevalence(df, threshold):

        if threshold == 0:
            print("Prevalence threshold was set to 0, thus no filtering will be applied.")
            return df

        # Pecentage of nonzero abundance of each genome
        float_cols       = df.select_dtypes(include='float').columns
        nonzero_fraction = (df[float_cols] != 0).sum(axis=1) / len(float_cols)

        # Split kept and removed rows
        df_filtered = df[nonzero_fraction >= threshold]
        df_removed  = df[nonzero_fraction < threshold]

        # Sum of filtered species kept
        sum_row = df_removed[float_cols].sum()
        new_row = {col: 0 for col in df.columns}

        new_row.update(sum_row.to_dict())

        new_row['classification'] = 'Removed_sum'
        new_row['genome']         = 'Total_removed_abd'

        # Append to filtered DataFrame
        df_filtered = pd.concat([df_filtered, pd.DataFrame([new_row])], ignore_index=True)

        return df_filtered


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Split dataset based on prevalence threshold.",
        usage="""
            ./split.py -a data/genome_counts_macro_filt_30_zerosrem.csv -g data/genome_metadata.csv -m data/sample_metadata_macro.csv \
                -b data/original_mets_data_20250113.xlsx -c day,treatment -p 0 -o microbetag_input/prev20 -s dpi,day_code,treatment_expl -n log1p
        """
    )
    parser.add_argument("-a", "--abundance-file", type=str, required=True, help="Path to the abundance file.")
    parser.add_argument("-g", "--genome-information", type=str, required=True, help="Path to the genome info file.")
    parser.add_argument("-b", "--metabolites-file", type=str, required=True, help="Path to the EXCEL file with metabolites.")
    parser.add_argument("-c", "--categories", type=lambda s: s.split(","), required=True, help="List of categories to consider.")
    parser.add_argument("-m", "--metadata-file", type=str, required=True, help="Path to the metadata file.")
    parser.add_argument("-o", "--outdir", type=str, required=False, default=None, help="Path to save output files.")
    parser.add_argument("-p", "--prevalence-threshold", type=float, required=False, default=0.2, help="Prevalence threshold (between 0 and 1).")
    parser.add_argument("-s", "--metadata-to-skip", type=lambda s: s.split(","), required=False, help="List of metadata columns to be removed from the metadata files to be built.")
    parser.add_argument("-n", "--normalize-metabolites", type=str, required=False, help="Normalization method for metabolics data. By default: None")

    args = parser.parse_args()

    splitter = SplitDataset(
        abundance_file   = args.abundance_file,
        genome_info      = args.genome_information,
        metadata_file    = args.metadata_file,
        metabolites_file = args.metabolites_file,
        categories       = args.categories,
        threshold        = args.prevalence_threshold,
        norm_metabo      = args.normalize_metabolites,
        outdir           = args.outdir,
        skip             = args.metadata_to_skip
    )

    splitter.overall()

    # splitter.split()
