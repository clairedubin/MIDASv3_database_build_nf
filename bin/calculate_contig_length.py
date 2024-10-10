#!/usr/bin/env python3
import argparse
from midas.common.utilities import get_contig_length
from midas.common.utils import OutputStream

def write_contig_length(dict_of_contig_length, contig_length_fp):
    with OutputStream(contig_length_fp) as stream:
        stream.write("\t".join(["genome_id", "contig_id", "contig_length"]) + "\n")
        for gid, r in dict_of_contig_length.items():
            for cid, clen in r.items():
                vals = [gid, cid, str(clen)]
                stream.write("\t".join(vals) + "\n")



if __name__ == "__main__":

    CLI=argparse.ArgumentParser()
    CLI.add_argument(
    "--fna_list",  
    nargs="*", 
    type=str,
    required=True,  
    )

    args = CLI.parse_args()

    fna_list = args.fna_list

    # We need this for within-contig-gene-position
    dict_of_contig_length = {}

    for genomefna in fna_list:
        genome_id = genomefna.replace('.fna','')
        dict_of_contig_length[genome_id] = get_contig_length(genomefna)

    with open("contigs.len", 'w') as f:
            f.write("\t".join(["genome_id", "contig_id", "contig_length"]) + "\n")
            for gid, r in dict_of_contig_length.items():
                for cid, clen in r.items():
                    vals = [gid, cid, str(clen)]
                    f.write("\t".join(vals) + "\n")
