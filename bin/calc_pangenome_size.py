#!/usr/bin/env python3
import argparse
import pandas as pd

if __name__ == "__main__":

    CLI=argparse.ArgumentParser()
    CLI.add_argument(
    "--gene_info_file",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--species",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--num_genomes",
    type=int,  
    required=True,
    )
    CLI.add_argument(
    "--rep_number",
    type=int,  
    required=True,
    )

    args = CLI.parse_args()
    gene_info_file = args.gene_info_file
    species = args.species
    num_genomes = args.num_genomes
    rep_number = args.rep_number


    rep_info = pd.Series({'species':species, 'num_genomes':num_genomes, 'rep_number':rep_number})

    df = pd.read_csv(gene_info_file, sep='\t')
    centroid_cols = [c for c in df.columns if 'centroid_' in c]
    counts = df.nunique().loc[centroid_cols].sort_index()

    merged = pd.concat([rep_info, counts]).to_frame().T
    merged.to_csv(f'pangenome_size_{species}_{num_genomes}_{rep_number}.csv', index=False)

