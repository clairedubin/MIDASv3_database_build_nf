#!/usr/bin/env python3
import os
from collections import defaultdict
import argparse
import pandas as pd
import re
from parse_centroids import read_uclust_info
from common import select_from_tsv

"""
Input:
    - temp/cdhit/gene_info.txt
    - temp/cdhit/centroids.99.ffn
    - temp/cdhit/genes.ffn
    - temp/cdhit/genes.len
Output:
    - genes_info.tsv: the unique id of this table is gene id of the cleaned set of genes
    - centroids.99.ffn
    - temp/centroids.xx.ffn
"""

GENE_LENGTH_SCHEMA = {
    "gene_id": str,
    "genome_id": str,
    "gene_length": int,
}

MARKER_INFO_SCHEMA = {
    "species_id": str,
    "genome_id": str,
    "gene_id": str,
    "gene_length": int,
    "marker_id": str
}


def scan_gene_length(gene_length_file):
    gene_length_dict = {}
    with open(gene_length_file) as stream:
        for r in select_from_tsv(stream, schema=GENE_LENGTH_SCHEMA, result_structure=dict):
            gene_length_dict[r["gene_id"]] = r["gene_length"]
    return gene_length_dict


def scan_mapfile(mapfile):
    """ Extract <marker_id, gene_id> pairs from the marker mapfile """
    dict_of_markers = {}
    with open(mapfile) as stream:
        # Example: 100001	GUT_GENOME000001	GUT_GENOME000001_01635	906	B000032
        for gene_id, marker_id in select_from_tsv(stream, ["gene_id", "marker_id"], schema=MARKER_INFO_SCHEMA):
            dict_of_markers[gene_id] = marker_id
    return dict_of_markers


def read_gene_info(centroid_info, gene_info_file, percent_id):
    # Parse intermediate gene_info.txt
    with open(gene_info_file) as stream:
        for r in select_from_tsv(stream, selected_columns=['gene_id', f'centroid_{percent_id}'], schema={'gene_id':str, f'centroid_{percent_id}':str}):
            centroid_info[r[0]][percent_id] = r[1]

def augment_gene_info(centroid_info, gene_to_marker, dict_of_gene_length, output_file_name):
    """ Augment gene_info.txt with two additional columns: gene_length and marker_id """

    df = pd.DataFrame.from_dict(centroid_info).T
    df.columns = [f'centroid_{i}' for i in df.columns]
    df.index.name = 'gene_id'
    df['gene_length'] = df.index.map(dict_of_gene_length)
    df['marker_id'] = df.index.map(gene_to_marker)
    df.to_csv(output_file_name, sep='\t')


def xref(cluster_files):
    # Again, let centroid_info[gene][percent_id] be the centroid of the percent_id
    # cluster containing gene; then reclustered to lower percent_id's.
    # Output: gene_info.txt for filtered centroids_99 genereated by cd-hit

    percents = cluster_files.keys()
    max_percent_id = max(percents)
    centroid_info = defaultdict(dict)
    for percent_id, uclust_file in cluster_files.items():
        if percent_id == max_percent_id:
            read_gene_info(centroid_info, uclust_file, percent_id)
        else:
            read_uclust_info(centroid_info, uclust_file, percent_id)

    # Check for a problem that occurs with improper import of genomes (when contig names clash).
    for g in centroid_info:

        cg = centroid_info[g][max_percent_id]
        ccg = centroid_info[cg][max_percent_id]

        assert cg == ccg, f"The {max_percent_id}-centroid relation should be idempotent, however, {cg} != {ccg}. {g}"

    # Infer coarser clusters assignments for all genes by transitivity
    for gc in centroid_info.values():
        gc_recluster = centroid_info[gc[max_percent_id]]
        for percent_id in percents:
            gc[percent_id] = gc_recluster[percent_id]

    return centroid_info

def split(iterable, portion_size):
        assert portion_size > 0
        p = []
        for item in iterable:
            p.append(item)
            if len(p) == portion_size:
                yield p
                p = []
        if p:
            yield p

if __name__ == "__main__":

    CLI=argparse.ArgumentParser()
    CLI.add_argument(
    "--gene_info_file",
    # nargs=1,
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--max_percent",
    # nargs=1,
    type=int,  
    required=True,
    )
    CLI.add_argument(
    "--gene_length_file",
    # nargs=1,
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--uclust_files",  
    nargs="*", 
    type=str,
    required=True,  
    )

    args = CLI.parse_args()

    gene_info_file = args.gene_info_file
    max_percent = args.max_percent
    gene_length_fp = args.gene_length_file
    uclust_files = args.uclust_files

    cluster_files = {}
    cluster_files[max_percent] = gene_info_file

    for f in sorted(uclust_files)[::-1]:

        assert os.path.exists(f)
        assert re.match(r'^uclust\.\d{1,2}\.txt$', f), \
        "File name does not match the required format uclust.XX.txt"

        cluster_percent = int(f.split('.')[1].replace('.txt',''))
        if cluster_percent == max_percent:
            continue

        cluster_files[cluster_percent] = f
    
    centroid_info = xref(cluster_files)

    dict_of_gene_length = scan_gene_length(gene_length_fp) # <gene_id:gene_length>
    gene_to_marker = scan_mapfile("temp/mapfile") # <gene_id, marker_id>

    # generate genes_info.tsv
    augment_gene_info(centroid_info, gene_to_marker, dict_of_gene_length, "genes_info.tsv")

