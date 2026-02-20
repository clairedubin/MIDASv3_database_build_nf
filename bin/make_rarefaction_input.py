#!/usr/bin/env python3

import argparse
import pandas as pd
import os

def add_genes_ffn(genome, species, db_path):
    return f'{db_path}/gene_annotations/{species}/{genome}/{genome}.genes.ffn'

def add_genes_len(genome, species, db_path):
    return f'{db_path}/gene_annotations/{species}/{genome}/{genome}.genes.len'


if __name__ == "__main__":

    CLI=argparse.ArgumentParser()
    CLI.add_argument(
    "--genomes_tsv",
    type=str,  
    required=True,
    )

    CLI.add_argument(
    "--species",
    type=str,  
    required=True,
    )

    CLI.add_argument(
    "--num_reps",
    type=int,  
    required=True,
    )
    
    CLI.add_argument(
    "--db_path",
    type=str,  
    required=True,
    )

    CLI.add_argument(
    "--max_size",
    type=int,
    required=False,
    help="Optional maximum size for sampling"
    )

    args = CLI.parse_args()

    genomes = args.genomes_tsv
    species = args.species
    num_reps = args.num_reps
    db_path = args.db_path

    df = pd.read_csv(genomes, sep='\t')
    df = df[df['species'].astype(str) == species]
    df = df[['genome','species']]

    df['genes_ffn'] = df.apply(lambda x: add_genes_ffn(x['genome'], x['species'], db_path), axis=1)
    df['genes_len'] = df.apply(lambda x: add_genes_len(x['genome'], x['species'], db_path), axis=1)
    df = df[(df['genes_ffn'].apply(os.path.exists)) & (df['genes_len'].apply(os.path.exists))]
    
    actual_rows = df.shape[0]

    if args.max_size is not None:
        max_size = args.max_size
    else:
        max_size = actual_rows

    start_val = max_size 
    
    i = start_val
    size_list = [start_val]

    while i >= 2:
       i *= 0.9
       size_list += [round(i)]

    for s in size_list:

        if s > actual_rows:
            print(f"Warning: Requested sample size {s} is larger than available genomes ({actual_rows}). Skipping.")
            continue

        if s == actual_rows:
            #we don't need to do this repeatedly since it is the max number of genomes and will yield the same result every time
            sample = df.sample(s, random_state=0)
            sample.to_csv(f'{species}_{s}genomes_rep0.csv', index=False)

        else:
            for rep in range(num_reps):
                sample = df.sample(s, random_state=s+rep)
                sample.to_csv(f'{species}_{s}genomes_rep{rep}.csv', index=False)
