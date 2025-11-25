#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import warnings
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
    "Adenine"
    "Guanosine"
    "Lactic acid"
    "Sulfoacetic acid"
]

class SplitDataset:
    def __init__(self, abundance_file, genome_info, metadata_file, metabolites_file, categories, threshold=0.2, outdir=None):

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
        metadata_f = metadata_f.drop(index=["batch", "tissue", "animal_m", "animal_y", "sample"], errors="ignore")

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

                # Pecentage of nonzero abundance of each genome
                float_cols       = tr_abd.select_dtypes(include='float').columns
                nonzero_fraction = (tr_abd[float_cols] != 0).sum(axis=1) / len(float_cols)

                # Split kept and removed rows
                df_filtered = tr_abd[nonzero_fraction >= self.threshold]
                df_removed  = tr_abd[nonzero_fraction < self.threshold]

                # sum of filtered species kept
                sum_row = df_removed[float_cols].sum()
                new_row = {col: 0 for col in tr_abd.columns}

                new_row.update(sum_row.to_dict())

                new_row['classification'] = 'Removed_sum'
                new_row['genome']         = 'Total_removed_abd'

                # Append to filtered DataFrame
                df_filtered = pd.concat([df_filtered, pd.DataFrame([new_row])], ignore_index=True)

                # Export those subset of the dataframes to files
                abd_outfile  = "_".join(["abd_prev", str(self.threshold), str(type), str(case)])
                abd_outfile += ".tsv"
                meta_outfile = "_".join(["metadata_prev", str(self.threshold), str(type), str(case)])
                meta_outfile += ".tsv"

                if self.outdir:
                    os.makedirs(self.outdir, exist_ok=True)
                df_filtered.to_csv(Path(self.outdir) / abd_outfile, index=False)
                tr_metadata.to_csv(Path(self.outdir) / meta_outfile, index_label="sample", sep="\t")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Split dataset based on prevalence threshold.")
    parser.add_argument("--abundance_file", type=str, required=True, help="Path to the abundance file.")
    parser.add_argument("--genome_information", type=str, required=True, help="Path to the genome info file.")
    parser.add_argument("--metabolites_file", type=str, required=True, help="Path to the EXCEL file with metabolites.")
    parser.add_argument("--categories", type=lambda s: s.split(","), required=True, help="List of categories to consider.")
    parser.add_argument("--metadata_file", type=str, required=True, help="Path to the metadata file.")
    parser.add_argument("--outdir", type=str, required=False, default=None, help="Path to save output files.")
    parser.add_argument("--prevalence_threshold", type=float, required=False, default=0.2, help="Prevalence threshold (between 0 and 1).")

    args = parser.parse_args()

    splitter = SplitDataset(
        abundance_file   = args.abundance_file,
        genome_info      = args.genome_information,
        metadata_file    = args.metadata_file,
        metabolites_file = args.metabolites_file,
        categories       = args.categories,
        threshold        = args.prevalence_threshold,
        outdir           = args.outdir
    )

    splitter.split()
