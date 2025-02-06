#!/usr/bin/env python3
import Bio.SeqIO
from common import has_ambiguous_bases
import argparse

if __name__ == "__main__":

    CLI=argparse.ArgumentParser()
    CLI.add_argument(
    "--genome",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--ffn",
    type=str,  
    required=True,
    )

    args = CLI.parse_args()
    ffn_path = args.ffn
    genome = args.genome

    # genome_id = "${genome}"
    output_genes = f"{genome}.genes.ffn"
    output_len = f"{genome}.genes.len"
    with open(output_genes, 'w') as o_genes, \
            open(output_len, 'w') as o_info:        
        
        for rec in Bio.SeqIO.parse(ffn_path, 'fasta'):
            gene_id = rec.id
            gene_seq = str(rec.seq).upper()
            gene_len = len(gene_seq)
            if gene_len <= 200 or has_ambiguous_bases(gene_seq) or gene_id == '' or gene_id == '|':
                pass
            else:
                o_genes.write(f">{gene_id}\n{gene_seq}\n")
                o_info.write(f"{gene_id}\t{genome}\t{gene_len}\n")