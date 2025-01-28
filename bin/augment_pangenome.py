#!/usr/bin/env python3
import os
import pandas as pd
import argparse
from functools import reduce

"""
INPUT: FROM build_pangenome + recluster_centroids
    - genes_info.tsv, produced by recluster_centroids.py
    - centroids.ffn, symlink to temp/cdhit/centroids.99.ffn
OUTPUT:
    - augment/{centroid_xx}_info.tsv: centroid_xx, centroid_xx_length, cluster_prevalence, cluster_gene_counts, cluster_marker_id, etc
    - augment/{centroid_xx}_lifted_to_{rep_genome}.tsv
    - augment/{centroid_xx}_stacked_matrix.tsv
"""

def list_cxx_coordinate_to_genome(genes_info, centroid_xx="centroid_99", qry_genome = ""):
    """ This function only locate the coordinate of centroids back to the query genome. """
    df = genes_info[['gene_id', centroid_xx]].copy()
    df['genome_id'] = df['gene_id'].apply(lambda x: extract_genomeid(x))
    cxx_lifted_to_qry = df[df['genome_id'] == qry_genome]
    cxx_lifted_to_qry = cxx_lifted_to_qry[['gene_id', centroid_xx]]
    # TODO: to identify multi-copy genes, we need to blast the lifetd pan-genes against qry genome
    return cxx_lifted_to_qry

def render_full_stacked_cxx_by_genome(genes_info, c_xx='centroid_99'):
    df = genes_info[['gene_id', c_xx]].copy()
    df['genome_id'] = df['gene_id'].apply(lambda x: extract_genomeid(x))
    centroid_xx_info = df[[c_xx, 'genome_id']].drop_duplicates()
    return centroid_xx_info

def extract_genomeid(s):
    # Replace everything after the right most "_" with empty string.
    index = s.rfind('_')
    return s[:index] if index != -1 else s

def generate_cluster_xx_info(genes_info, c_xx = 'centroid_99', cutoff=0.5):
    # Say we start with c99, which is the original cluster_info
    c_xx_df = genes_info[genes_info["gene_id"] == genes_info[c_xx]]
    c_xx_prev = compute_cxx_prevalence(genes_info, c_xx)
    c_xx_gene_counts = compute_cxx_gene_counts(genes_info, c_xx)
    c_xx_length = compute_cxx_length(genes_info, c_xx, 'max') #<---
    c_xx_markers = impute_cxx_marker_id(genes_info, c_xx, cutoff)
    #c_xx_gene_info = decorate_cxx_with_gene_info(genes_info, c_xx)
    if c_xx == 'centroid_99':
        c_xx_df = c_xx_df.drop(['gene_id', 'gene_length', 'marker_id'], axis=1)
    else:
        c_xx_df = c_xx_df[c_xx]
    list_of_dfs = [c_xx_df, c_xx_prev, c_xx_gene_counts, c_xx_length, c_xx_markers]
    result = reduce(lambda left, right: pd.merge(left, right, on=c_xx, how='left'), list_of_dfs)
    # Replace NaN with 0 for marker_density
    result[f'{c_xx}_marker_density'] = result[f'{c_xx}_marker_density'].fillna(0)
    return result

def compute_cxx_gene_counts(genes_info, centroid_xx="centroid_99"):
    c_xx_gene_counts = genes_info[[centroid_xx, 'gene_id']].groupby(centroid_xx).size().reset_index(name='total_gene_counts')
    c_xx_gene_counts.columns = [centroid_xx, f'{centroid_xx}_gene_counts']
    return c_xx_gene_counts

def compute_cxx_length(genes_info, centroid_xx = "centroid_99", func='median'):
    """ For each centroid_xx, we compute either (1) median gene length or (2) max gene length across all gene members """
    if func == "median":
        cxx_length = genes_info[[centroid_xx, 'gene_length']].groupby(centroid_xx)['gene_length'].median().reset_index()
    else:
        cxx_length = genes_info[[centroid_xx, 'gene_length']].groupby(centroid_xx)['gene_length'].max().reset_index()
    cxx_length.columns = [centroid_xx, f'{centroid_xx}_gene_length']
    return cxx_length

def compute_cxx_prevalence(genes_info, c_xx='centroid_99'):
    """ For each centroid_xx, we compute the operational gene family / cluster prevalence across all genomes """
    # Parse genome_id from the gene_id
    df = genes_info[['gene_id', c_xx]].copy()
    df['genome_id'] = df['gene_id'].apply(lambda x: extract_genomeid(x))
    # Subset the genes_info with custom list of genomes
    total_genomes = df['genome_id'].nunique()
    # Group_by centroid_xx and count the number of genomes with genes falling into given operational gene cluster
    centroid_xx_info = df[[c_xx, 'genome_id']].drop_duplicates()
    centroid_xx_prev = centroid_xx_info.groupby(c_xx)['genome_id'].count().reset_index(name='genome_counts')
    centroid_xx_prev['prevalence'] = centroid_xx_prev['genome_counts'] / total_genomes
    centroid_xx_prev.columns = [c_xx, f"{c_xx}_genome_counts", f"{c_xx}_genome_prevalence"]
    return centroid_xx_prev

def impute_cxx_marker_id(genes_info, c_xx="centroid_99", cutoff=0.5):
    """ Marker assignment for each centroid_xx is voted by all gene members """
    df = genes_info[[c_xx, 'marker_id']]
    # First compute the total_gene_counts per centroid_xx
    gene_counts = df.groupby(c_xx).size().reset_index(name='total_gene_counts')
    # Second compute the centroid_xx, marker_id occurrences, exlucding NA marker_ids
    pair_counts = df.dropna(subset=['marker_id']).groupby([c_xx, 'marker_id']).size().reset_index(name='occurrence')
    max_pair_counts = select_marker_by_max_counts(pair_counts, c_xx)
    # Third compute the ratio
    c_xx_markers = pd.merge(max_pair_counts, gene_counts, on=c_xx)
    c_xx_markers['marker_density'] = c_xx_markers['occurrence'] / c_xx_markers['total_gene_counts']
    # Fouth only keep marker assignment if (max) occurrence > cutoff
    passed_markers = c_xx_markers[c_xx_markers['marker_density'] > cutoff][[c_xx, "marker_id"]]
    # Fifth we still want to keep the low marker_density
    c_xx_markers = pd.merge(c_xx_markers[[c_xx, 'marker_density']], passed_markers, on=c_xx, how='left')
    c_xx_markers.columns = [c_xx, f'{c_xx}_marker_density', f'{c_xx}_marker_id']
    return c_xx_markers

def select_marker_by_max_counts(df, centroid_xx):
    # This is more like a remedy than a solution.
    # Edge case 1: cxx has two both marker_id and empty marker_id (NaN)
    # Edge case 2: rarely a erroneous fusion ORF could be assigned to two markers, and we keep the marker_id assignment with max gene_counts
    # Edge case 3: we drop ties
    df['max_counts'] = df.groupby(centroid_xx)['occurrence'].transform('max')
    max_occur_df = df[df['occurrence'] == df['max_counts']]
    max_occur_df = max_occur_df.drop(columns=['max_counts'])
    max_occur_df = max_occur_df[~max_occur_df.duplicated(subset=centroid_xx, keep=False)]
    return max_occur_df

def write_cxx_info(c99_df, cluster_info_fp):
    c99_df.to_csv(cluster_info_fp, sep='\t', index=False)

if __name__ == "__main__":

    CLI=argparse.ArgumentParser()
    CLI.add_argument(
    "--gene_info_file",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--cluster_thresholds",
    nargs="*",
    type=int,  
    required=True,
    )
    CLI.add_argument(
    "--rep_genome",
    type=str,  
    required=True,
    )

    args = CLI.parse_args()

    genes_info_fp = args.gene_info_file
    cluster_thresholds = args.cluster_thresholds
    rep_genome= args.rep_genome

    assert os.path.isfile(genes_info_fp), f"Pangenome output {genes_info_fp} are missing."
    genes_info = pd.read_csv(genes_info_fp, delimiter="\t")

    for xx in cluster_thresholds:
        centroid_xx = f"centroid_{xx}"
        cluster_xx_info = generate_cluster_xx_info(genes_info, centroid_xx)
        write_cxx_info(cluster_xx_info, f"clusters_{xx}_info.tsv")

    lowest_threshold, highest_threshold = min(cluster_thresholds), max(cluster_thresholds)

    centroid_xx_info = render_full_stacked_cxx_by_genome(genes_info, f"centroid_{lowest_threshold}")
    centroid_xx_info.to_csv(f"centroid_{lowest_threshold}_stacked_matrix.tsv", sep='\t', index=False)

    cxx_lifted_to_qry = list_cxx_coordinate_to_genome(genes_info, f"centroid_{highest_threshold}", rep_genome)
    cxx_lifted_to_qry.to_csv(f"centroid_{highest_threshold}_lifted_to_{rep_genome}.tsv", sep='\t', index=False)