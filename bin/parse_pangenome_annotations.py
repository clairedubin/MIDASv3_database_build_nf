#!/usr/bin/env python3
import os
import sys
from itertools import chain
import pandas as pd
import argparse
from pybedtools import BedTool

from midas.common.argparser import add_subcommand
from midas.common.utils import tsprint, command, multithreading_map, pythonpath, num_physical_cores, copy_star
from midas.common.utilities import decode_genomes_arg, decode_species_arg, scan_eggnog
from midas.models.midasdb import MIDAS_DB
from midas.params.inputs import MIDASDB_NAMES

"""
For given genome, we take the following inputs:
Input:
    - Species: pangenomes/cluster_info.txt, contigs_len
    - Genome: genes.feature
    - Funannot: genomad, mefinder and resfinder results
Output:
    - Four headerless temp tables for centroids.99
    """
def parse_genomad_virus_genes(file_path):
    """ Parse Genomad virus genes TSV into DataFrame """
    df = pd.read_csv(file_path, delimiter="\t")

    if df.empty:
        return pd.DataFrame()

    df = df[['gene', 'start', 'end', 'annotation_conjscan', 'annotation_amr', 'annotation_accessions', 'annotation_description']]
    df[['contig_plus', 'gene_num']] = df['gene'].str.rsplit('_', n=1, expand=True)

    # For the provisus case, we only keep the contig_id
    if df['contig_plus'].str.contains(r'\|p').any():
        split_data = df['contig_plus'].str.split(r'\|p', expand=True)
        # If the delimiter is not present, then the second column will be NaN.
        # We'll use that to replace NaNs in the first column with the original 'contig_plus' values.
        split_data[0].where(split_data[1].notna(), df['contig_plus'], inplace=True)
        # Assign the results back to the original dataframe
        df['contig_external'] = split_data[0]
    else:
        df['contig_external'] = df['contig_plus']

    assert ~df[['contig_external', 'start', 'end']].duplicated().any(), f"Duplicated virus results for {file_path}"

    df = df[['contig_external', 'start', 'end', 'gene',  'annotation_conjscan', 'annotation_amr', 'annotation_accessions', 'annotation_description']]
    df = df.rename(columns={'start': 'start_genomad', 'end': 'end_genomad', 'gene':'gene_genomad'})

    return df


def parse_genomad_plasmid_genes(file_path):
    """ Parse Genomad plasmid genes TSV into DataFrame """
    df = pd.read_csv(file_path, delimiter="\t")

    if df.empty:
        return pd.DataFrame()

    df = df[['gene', 'start', 'end', 'annotation_conjscan', 'annotation_amr', 'annotation_accessions', 'annotation_description']]
    df[['contig_external', 'gene_num']] = df['gene'].str.rsplit('_', n=1, expand=True)

    assert ~df[['contig_external', 'start', 'end']].duplicated().any(), f"Duplicated plasmid results for {file_path}"

    df = df[['contig_external', 'start', 'end', 'gene', 'annotation_conjscan', 'annotation_amr', 'annotation_accessions', 'annotation_description']]
    df = df.rename(columns={'start': 'start_genomad', 'end': 'end_genomad', 'gene': 'gene_genomad'})

    return df


def parse_mefinder(file_path):
    """ Parse MEFINDER into DataFrame """
    df = pd.read_csv(file_path, delimiter=",", skiprows=5)

    if df.empty:
        return pd.DataFrame()

    df['contig_id'] = df['contig'].str.split(' ').str[0]
    df = df[['contig_id', 'start', 'end', 'mge_no', 'prediction', 'name', 'type', 'synonyms']]
    df = df.rename(columns={'contig_id': 'contig_external', 'start': 'start_mefinder', 'end': 'end_mefinder'})

    return df


def parse_resfinder(file_path):
    # rename the columns without spaces
    new_columns = ['resistance_gene', 'identity', 'align_ratio', 'coverage', 'within_reference_position', 'contig', 'within_contig_position', 'phenotype', 'accession_no']
    df = pd.read_csv(file_path, sep='\t', header=0, names=new_columns)

    if df.empty:
        return pd.DataFrame()

    df['contig'] = df['contig'].str.split(' ').str[0]
    df[['start', 'end']]  = df['within_contig_position'].str.split('\\.\\.', expand=True)
    # Notes there can be duplicated annotaitons for the same gene range (contig-start-end)
    df = df[['contig', 'start', 'end', 'resistance_gene', 'phenotype', 'accession_no']]
    df = df.rename(columns={'contig': 'contig_external', 'start': 'start_resfinder', 'end': 'end_resfinder'})
    return df


def merge_annot_with_genes(df, all_genes):
    bed2 = BedTool.from_dataframe(df)
    bed1 = BedTool.from_dataframe(all_genes)
    overlaps = bed1.intersect(bed2, wa=True, wb=True)

    if len(overlaps) == 0:
        return pd.DataFrame()

    overlapping_df = pd.read_table(overlaps.fn, header=None)
    overlapping_df.columns = list(all_genes.columns) + list(df.columns)
    # Drop duplicated column: contig_id
    overlapping_df = overlapping_df.drop(columns=['contig_external'])
    # Reorder the column names
    overlapping_df = overlapping_df[['gene_id'] + [col for col in overlapping_df if col != 'gene_id']]
    return overlapping_df


def write_processed_genome_annot(df, all_genes, local_dest):
    if df.empty:
        command(f"touch {local_dest}")
    else:
        df = merge_annot_with_genes(df, all_genes)
        if df.empty:
            command(f"touch {local_dest}")
        else:
            df.to_csv(local_dest, sep='\t', index=False, header=False)

if __name__ == "__main__":

    CLI=argparse.ArgumentParser()
    CLI.add_argument(
    "--species",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--genome",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--max_cluster_info_file",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--contig_len_file",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--genes_file",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--genomad_virus_file",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--genomad_plasmid_file",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--mefinder_csv",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--resfinder_txt",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--eggnog_results_file",
    type=str,  
    required=True,
    )

    args = CLI.parse_args()

genome_id = args.genome
species_id = args.species
cluster_info_fp = args.max_cluster_info_file
contig_len_fp = args.contig_len_file
gene_feature_fp = args.genes_file

contig_len = pd.read_csv(contig_len_fp, sep='\t')
gene_features = pd.read_csv(gene_feature_fp, sep='\t')

# 2023-12-05: we keep ALL the GENES, instead of centroids_99
all_genes = pd.merge(gene_features, contig_len[['contig_id', 'contig_length']], left_on='contig_id', right_on='contig_id', how='inner')
all_genes = all_genes[['contig_id', 'start', 'end', 'gene_id', 'strand', 'gene_type', 'contig_length']]

# Genomad virus
df = parse_genomad_virus_genes(args.genomad_virus_file)
local_dest = f"genomad_virus/{genome_id}_genomad_virus.tsv"
write_processed_genome_annot(df, all_genes, local_dest)

# Genomad plasmid
df = parse_genomad_plasmid_genes(args.genomad_plasmid_file)
local_dest = f"genomad_plasmid/{genome_id}_genomad_plasmid.tsv"
write_processed_genome_annot(df, all_genes, local_dest)

# MEfinder
df = parse_mefinder(args.mefinder_csv)
local_dest = f"mefinder/{genome_id}_mefinder.tsv"
write_processed_genome_annot(df, all_genes, local_dest)

# Resfinder
df = parse_resfinder(args.resfinder_txt)
local_dest = f"resfinder/{genome_id}_resfinder.tsv"
write_processed_genome_annot(df, all_genes, local_dest)

# TODO: make it so centroids_99 is not hard coded below, eg. centroids_{max_percent}

# # EggNOG was run on the centroid_99s.
centroids_99 = pd.read_csv(cluster_info_fp, sep='\t', usecols=[0]) # Only keep the list of centroid_99
centroids_99 = pd.merge(centroids_99, gene_features, left_on='centroid_99', right_on='gene_id', how='inner')
centroids_99 = pd.merge(centroids_99, contig_len[['contig_id', 'contig_length']], left_on='contig_id', right_on='contig_id', how='inner')
centroids_99 = centroids_99[['contig_id', 'start', 'end', 'centroid_99', 'strand', 'gene_type', 'contig_length']]

df = scan_eggnog(args.eggnog_results_file)
local_dest = f"eggnog/{genome_id}_eggnog.tsv"
merged_df  = df.merge(centroids_99, left_on='#query', right_on='centroid_99', how='inner')
merged_df = merged_df.drop(columns=['centroid_99'])
merged_df.to_csv(local_dest, sep='\t', index=False, header=False)
