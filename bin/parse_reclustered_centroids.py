#!/usr/bin/env python3
import os
from collections import defaultdict
import argparse
import pandas as pd
from midas.common.utils import select_from_tsv, InputStream, OutputStream, cat_files
from midas.params.schemas import GENE_LENGTH_SCHEMA, MARKER_INFO_SCHEMA, PANGENOME_INFO_SCHEMA


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
###TODO: make these shared functions with parse_centroids.py

def parse_uclust(uclust_file, select_columns):
    # The uclust TSV file does not contain a header line.  So, we have to hardcode the schema here.  Then select specified columns.
    all_uclust_columns = ['type', 'cluster_id', 'size', 'pid', 'strand', 'skip1', 'skip2', 'skip3', 'gene_id', 'centroid_id']
    with open(uclust_file, 'r') as ucf:
        for r in select_from_tsv(ucf, select_columns, all_uclust_columns):
            yield r


def read_uclust_info(centroid_info, uclust_file, percent_id):
    # Get centroid_info from uclust
    for r_type, r_gene, r_centroid in parse_uclust(uclust_file, ['type', 'gene_id', 'centroid_id']):
        if r_type == 'S':
            # r itself is the centroid of its cluster
            centroid_info[r_gene][percent_id] = r_gene
        elif r_type == 'H':
            # r is not itself a centroid
            centroid_info[r_gene][percent_id] = r_centroid
        else:
            # ignore all other r types
            pass


# @retry
def scan_gene_length(gene_length_file):
    gene_length_dict = {}
    with InputStream(gene_length_file) as stream:
        for r in select_from_tsv(stream, schema=GENE_LENGTH_SCHEMA, result_structure=dict):
            gene_length_dict[r["gene_id"]] = r["gene_length"]
    return gene_length_dict


# @retry
def scan_mapfile(mapfile):
    """ Extract <marker_id, gene_id> pairs from the marker mapfile """
    dict_of_markers = {}
    with InputStream(mapfile) as stream:
        # Example: 100001	GUT_GENOME000001	GUT_GENOME000001_01635	906	B000032
        for gene_id, marker_id in select_from_tsv(stream, ["gene_id", "marker_id"], schema=MARKER_INFO_SCHEMA):
            dict_of_markers[gene_id] = marker_id
    return dict_of_markers


def read_gene_info(centroid_info, gene_info_file, percent_id):
    # Parse intermediate gene_info_cdhit.tsv
    with InputStream(gene_info_file) as stream:
        for r in select_from_tsv(stream, selected_columns=['gene_id', f'centroid_{percent_id}'], schema={'gene_id':str, f'centroid_{percent_id}':str}):
            centroid_info[r[0]][percent_id] = r[1]


def augment_gene_info(centroid_info, gene_to_marker, dict_of_gene_length, output_file_name):
    """ Augment gene_info.txt with two additional columns: gene_length and marker_id """

    df = pd.DataFrame.from_dict(centroid_info).T
    df['gene_length'] = df.index.map(dict_of_gene_length)
    df['marker_id'] = df.index.map(gene_to_marker)
    df.to_csv(output_file_name)

    ## redone as above because if centroid cols weren't in order, they would assign the wrong labels
    # with OutputStream(output_file_name) as stream:
    #     stream.write("\t".join(['gene_id']+centroid_cols+["gene_length","marker_id"]) + "\n")
    #     for gene_id, r in centroid_info.items():
    #         gene_len = dict_of_gene_length[gene_id]
    #         marker_id = gene_to_marker[gene_id] if gene_id in gene_to_marker else ""
    #         val = [gene_id] + list(r.values()) + [gene_len, marker_id]
    #         stream.write("\t".join(map(str, val)) + "\n")


def xref(cluster_files):
    # Again, let centroid_info[gene][percent_id] be the centroid of the percent_id
    # cluster containing gene; then reclustered to lower percent_id's.
    # Output: gene_info.txt for filtered centroids_99 genereated by cd-hit

    percents = cluster_files.keys()
    max_percent_id = max(percents)
    centroid_info = defaultdict(dict)
    for percent_id, (_, uclust_file) in cluster_files.items():
        if percent_id == max_percent_id:
            read_gene_info(centroid_info, uclust_file, percent_id)
        else:
            read_uclust_info(centroid_info, uclust_file, percent_id)

    # Check for a problem that occurs with improper import of genomes (when contig names clash).
    for g in centroid_info:
        cg = centroid_info[g][max_percent_id]
        ccg = centroid_info[cg][max_percent_id]
        assert cg == ccg, f"The {max_percent_id}-centroid relation should be idempotent, however, {cg} != {ccg}."

    # Infer coarser clusters assignments for all genes by transitivity
    for gc in centroid_info.values():
        gc_recluster = centroid_info[gc[max_percent_id]]
        for percent_id in percents:
            gc[percent_id] = gc_recluster[percent_id]

    return centroid_info


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
    CLI.add_argument(
    "--genome_marker_files",
    nargs="*",
    type=str,  
    required=True,
    )

    # parse the command line
    args = CLI.parse_args()

    gene_info_file = args.gene_info_file
    max_percent = args.max_percent
    gene_length_fp = args.gene_length_file
    marker_map_files = args.genome_marker_files
    uclust_files = args.uclust_files

    cluster_files = {}
    cluster_files[max_percent] = ['', gene_info_file]

    for f in sorted(uclust_files):
        assert os.path.exists(f)
        ##TODO: assert that f matches a regex uclust.XX.txt

        cluster_percent = int(f.split('.')[1].replace('.txt',''))

        ##including empty str so dict format is compatible with definitions above
        ##can change this eventually
        cluster_files[cluster_percent] = ['', f]

    centroid_info = xref(cluster_files)
    print(centroid_info)

    # augment temp/vsearch/gene_info.txt with gene_length
    # gene_length_fp = midas_db.get_target_layout("pangenome_tempfile", False, species_id, "vsearch", "genes.len")
    dict_of_gene_length = scan_gene_length(gene_length_fp) # <gene_id:gene_length>

    # collect SGC mapfile from all genomes of species_id
    mapfiles_by_genomes = {}
    for f in marker_map_files:
        assert os.path.exists(f)
        genome_name = f.replace('.markers.map','')
        mapfiles_by_genomes[genome_name] = f

    cat_files(mapfiles_by_genomes.values(), "mapfile", 20)
    gene_to_marker = scan_mapfile("mapfile") # <gene_id, marker_id>

    # generate genes_info.tsv
    augment_gene_info(centroid_info, gene_to_marker, dict_of_gene_length, "genes_info.tsv")

