#!/usr/bin/env python3
import os
import sys
import argparse
from midas.common.argparser import add_subcommand, SUPPRESS
from midas.common.utils import tsprint, InputStream, retry, command, multithreading_map, find_files, upload, pythonpath, upload_star, num_physical_cores
from midas.common.utilities import scan_genes, decode_genomes_arg
from midas.models.midasdb import MIDAS_DB
from midas.params.inputs import hmmsearch_max_evalue, hmmsearch_min_cov, MIDASDB_NAMES

###EXAMPLE USAGE

#python3 bin/infer_markers.py --genome GCA_007120565.1 --species 117088 --hmmsearch_file test_db/markers/phyeco/temp/117088/GCA_007120565.1/GCA_007120565.1.hmmsearch --annotation_ffn test_db/gene_annotations/117088/GCA_007120565.1/GCA_007120565.1.ffn

def find_hits(hmmsearch_file):
    hits = {}
    for r in parse_hmmsearch(hmmsearch_file):
        if r['evalue'] > hmmsearch_max_evalue:
            continue
        if min(r['qcov'], r['tcov']) < hmmsearch_min_cov:
            continue
        if r['target'] not in hits:
            hits[r['target']] = r
        elif r['evalue'] < hits[r['target']]['evalue']:
            hits[r['target']] = r
    return list(hits.values())

def parse_hmmsearch(hmmsearch_file):
    """ Parse HMMER domblout files. Return data-type formatted dictionary """
    with InputStream(hmmsearch_file) as f_in:
        for line in f_in:
            if line[0] == "#":
                continue
            x = line.rstrip().split()
            query = x[0]
            target = x[3]
            evalue = float(x[12])
            qcov = (int(x[20]) - int(x[19]) + 1)/float(x[2])
            tcov = (int(x[16]) - int(x[15]) + 1)/float(x[5])
            yield {'query':query, 'target':target, 'evalue':evalue, 'qcov':qcov, 'tcov':tcov, 'qlen':int(x[2]), 'tlen':int(x[5])}

def compute_marker_genes(args):
    """ Parse local HMM search output file """
    hmmsearch_file = args.hmmsearch_file
    annotation_ffn = args.annotation_ffn
    
    genes = scan_genes(annotation_ffn)

    hmmsearch_seq = f"{args.genome}.markers.fa"
    hmmsearch_map = f"{args.genome}.markers.map"

    with open(hmmsearch_seq, "w") as o_seq, open(hmmsearch_map, "w") as o_map:
        for rec in find_hits(hmmsearch_file):
            marker_gene = genes[rec["query"]].upper()
            marker_info = [args.species, args.genome, rec["query"], len(marker_gene), rec["target"]]
            o_map.write('\t'.join(str(mi) for mi in marker_info) + '\n')
            o_seq.write('>%s\n%s\n' % (rec['query'], marker_gene))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome',
                           dest='genome',
                           required=True,
                           help="genome ID")
    parser.add_argument('--species',
                           dest='species',
                           required=True,
                           help="species ID")
    parser.add_argument('--hmmsearch_file',
                           dest='hmmsearch_file',
                           required=True,
                           help="path to hmmsearch_file") 
    parser.add_argument('--annotation_ffn',
                           dest='annotation_ffn',
                           required=True,
                           help="path to annotation_ffn") 
    args = parser.parse_args()
    compute_marker_genes(args)
