import argparse
import pandas as pd
import numpy as np
from io import StringIO

def decorate_genes_info_with_annot(genes_info, annot_files):
    # decorate genes_info with binary gene annotation
    for mstep, mfile in annot_files.items():
        new_col_name = fetch_new_col_name(mstep)
        annot_df = pd.read_csv(mfile, sep='\t')
        annot_df[new_col_name] = 1
        annot_df = annot_df[['gene_id', new_col_name]].groupby('gene_id').first().reset_index()
        # append gene level MGE status
        genes_info = pd.merge(genes_info, annot_df, on='gene_id', how='left')
        genes_info[new_col_name] = np.where(genes_info[new_col_name].isna(), 0, 1)
    return genes_info

def fetch_new_col_name(mstep, by="gene"):
    if mstep == "genomad_virus":
        new_col_name = f"{by}_is_phage"
    if mstep == "genomad_plasmid":
        new_col_name = f"{by}_is_plasmid"
    if mstep == "mefinder":
        new_col_name = f"{by}_is_amr"
    if mstep == "resfinder":
        new_col_name = f"{by}_is_me"
    return new_col_name

def scan_eggnog(filename):
    # Read the file line-by-line and exclude lines starting with '##'
    with open(filename, 'r') as file:
        lines = [line for line in file if not line.startswith('##')]
    # Convert the filtered lines into a DataFrame
    df = pd.read_csv(StringIO('\n'.join(lines)), sep='\t')
    return df


def annotation_ratio_x_members(genes_annotated, eggnog_file, xx='75'):
    """ Compuate the annotation ratio for centroid_xx level across its all gene members' annotation """
    eggnog = scan_eggnog(eggnog_file)
    eggnog.drop_duplicates(inplace=True)

    centroids_xx = f'centroid_{xx}'
    target_cols = ['gene_is_phage', 'gene_is_plasmid', 'gene_is_amr', 'gene_is_me']
    cxx_df = genes_annotated.groupby([centroids_xx])[target_cols].sum() / genes_annotated.groupby([centroids_xx])[target_cols].count()
    cxx_df.columns = ['phage_ratio', 'plasmid_ratio', 'amr_ratio', 'me_ratio']
    cxx_df = cxx_df.reset_index()

    eggnog_cols = ['#query', 'COG_category', 'EC', 'KEGG_Pathway', 'Description', 'PFAMs']
    cxx_df = pd.merge(cxx_df, eggnog[eggnog_cols], left_on = centroids_xx, right_on = '#query', how='left')
    cxx_df = cxx_df.drop(columns=['#query'])
    return cxx_df

if __name__ == "__main__":

    CLI=argparse.ArgumentParser()
    CLI.add_argument(
    "--cluster_thresholds",
    nargs="*",
    type=int,  
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
    "--mefinder_file",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--resfinder_file",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--eggnog_file",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--genes_info_file",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--cluster_info_file_list",
    nargs="*",
    type=str,  
    required=True,
    )
    
    
    args = CLI.parse_args()
    cluster_thresholds = args.cluster_thresholds

    genes_info_file = args.genes_info_file

    annotation_tools = ['genomad_virus', 'genomad_plasmid', 'mefinder', 'resfinder', 'eggnog']
    annot_files = {
        'genomad_virus':args.genomad_virus_file,
        'genomad_plasmid':args.genomad_plasmid_file, 
        'mefinder':args.mefinder_file,
        'resfinder':args.resfinder_file,
    }

    genes_info = pd.read_csv(genes_info_file, sep="\t")
    genes_annotated = decorate_genes_info_with_annot(genes_info, annot_files)
    genes_annotated.to_csv("genes_annotated.tsv", sep='\t', index=False)

    # Compute centroid_xx annotation density
    for xx in cluster_thresholds:
        cxx_df = annotation_ratio_x_members(genes_annotated, args.eggnog_file, xx)
        cxx_df.to_csv(f"clusters_{xx}_annot.tsv",  sep='\t', index=False)

    # Merge the cluster info with annotation
    for xx in cluster_thresholds:
        xx_info = pd.read_csv(f"old/clusters_{xx}_info.tsv", sep="\t")
        xx_anno = pd.read_csv(f"clusters_{xx}_annot.tsv", sep="\t")
        assert xx_info.shape[0] == xx_anno.shape[0], f"The shape of augment/clusters_{xx}_info.tsv disagree with annotation/clusters_{xx}_annoted.tsv"
        xx_df = pd.merge(xx_info, xx_anno, on=f"centroid_{xx}", how='left')
        xx_df.to_csv(f"clusters_{xx}_info.tsv", sep='\t', index=False)