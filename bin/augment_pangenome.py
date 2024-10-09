#!/usr/bin/env python3
import os
import pandas as pd
import argparse

from midas.common.utils import OutputStream
from midas.common.utilities import render_full_stacked_cxx_by_genome, generate_cluster_xx_info, list_cxx_coordinate_to_genome

"""
INPUT::FROM build_pangenome + recluster_centroids
    - genes_info.tsv, produced by recluster_centroids.py
    - centroids.ffn, symlink to temp/cdhit/centroids.99.ffn
OUTPUT:
    - augment/{centroid_xx}_info.tsv: centroid_xx, centroid_xx_length, cluster_prevalence, cluster_gene_counts, cluster_marker_id, etc
    - augment/{centroid_xx}_lifted_to_{rep_genome}.tsv
    - augment/{centroid_xx}_stacked_matrix.tsv
"""

def write_contig_length(dict_of_contig_length, contig_length_fp):
    with OutputStream(contig_length_fp) as stream:
        stream.write("\t".join(["genome_id", "contig_id", "contig_length"]) + "\n")
        for gid, r in dict_of_contig_length.items():
            for cid, clen in r.items():
                vals = [gid, cid, str(clen)]
                stream.write("\t".join(vals) + "\n")


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