#!/usr/bin/env python3
import gffutils
import argparse

if __name__ == "__main__":

    CLI=argparse.ArgumentParser()
    CLI.add_argument(
    "--gff",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--genome",
    type=str,  
    required=True,
    )

    args = CLI.parse_args()
    gff_path = args.gff
    genome = args.genome

    db = gffutils.create_db(gff_path, f'{gff_path}.db')
    to_write = "\t".join(["gene_id", "contig_id", "start", "end", "strand", "gene_type"]) + "\n"
    for feature in db.all_features():
        if feature.source == "prokka":
            continue
        if "ID" not in feature.attributes: 
            continue
        seqid = feature.seqid
        start = feature.start
        stop = feature.stop
        strand = feature.strand
        gene_id = feature.attributes['ID'][0]
        locus_tag = feature.attributes['locus_tag'][0]
        assert gene_id == locus_tag
        gene_type = feature.featuretype
        to_write += "\t".join([gene_id, seqid, str(start), str(stop), strand, gene_type]) + "\n"

    with open(f'{genome}.genes', 'w') as f:
        f.write(to_write)