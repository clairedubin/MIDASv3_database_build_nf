#!/usr/bin/env python3

from collections import defaultdict
from math import floor
import Bio
from Bio import SeqIO
import argparse
import json

def scan_fasta(fasta_file):
    """ Scan FASTA FILE to get seq and len """
    seqs = {}
    for rec in Bio.SeqIO.parse(fasta_file, 'fasta'):
        seqs[rec.id] = {
            "id": rec.id,
            "length": len(rec.seq),
            "seq": str(rec.seq),
        }
    return seqs

def design_run_snps_chunks(species_id, contigs_file, chunk_size):
    """ Given the Genome and chunk_size, the structure of the chunks are the same.
        Each chunk is indexed by (species_id, chunk_id) """

    # The structure of the chunks depends on the representative genome sequences
    contigs = scan_fasta(contigs_file)

    # Start with full chunks
    chunk_id = 0
    chunks_of_sites = defaultdict(list) # list of tuples
    unassigned_contigs = defaultdict(dict)
    max_contig_length = 0

    for contig in contigs.values():
        contig_length = contig["length"]
        if contig_length > max_contig_length:
            max_contig_length = contig_length

        contig_id = contig["id"]
        # left closed, right open
        if contig_length < chunk_size:
            unassigned_contigs[contig_id] = {"contig_id": contig_id,
                                             "contig_start": 0,
                                             "contig_end": contig_length,
                                             "contig_length": contig_length,
                                             "compute_reads": True}
        else:
            number_of_full_chunks = floor(contig_length/chunk_size)
            for ni, ci in enumerate(range(0, contig_length, chunk_size)):
                if ni == number_of_full_chunks: # last chunk
                    unassigned_contigs[contig_id] = {"contig_id": contig_id,
                                                     "contig_start": ci,
                                                     "contig_end": contig_length,
                                                     "contig_length": contig_length - ci,
                                                     "compute_reads": False}
                else:
                    count_flag = ni == 0 # first chunk
                    chunks_of_sites[chunk_id] = [(species_id, chunk_id, contig_id, ci, ci+chunk_size, count_flag, 0)]
                    chunk_id += 1

    if unassigned_contigs:
        # Partition unassigned short contigs into subsets
        subset_of_contigs, chunk_id = partition_contigs_into_chunks(unassigned_contigs, chunk_size, chunk_id)

        # Add the partitioned subsets to chunks
        for chunk_dict in subset_of_contigs.values():
            _chunk_id = chunk_dict["chunk_id"]
            list_of_contigs = chunk_dict["contigs_id"]
            for _cidx, _cid in enumerate(list_of_contigs):
                cstart = unassigned_contigs[_cid]["contig_start"]
                cend = unassigned_contigs[_cid]["contig_end"]
                cflag = unassigned_contigs[_cid]["compute_reads"]
                chunks_of_sites[_chunk_id].append((species_id, _chunk_id, _cid, cstart, cend, cflag, _cidx))
        assert chunk_id == _chunk_id+1

    # Finally the merge jobs
    number_of_chunks = chunk_id
    chunks_of_sites[-1] = (species_id, -1, number_of_chunks, max_contig_length)

    return chunks_of_sites

def partition_contigs_into_chunks(unassigned_contigs, chunk_size, chunk_id):
    """ Partition short, unassigned contigs into subsets/chunks.
        Similar to the problem of partition to K equal sum subsets """

    # Sort contigs by descending order of contig length
    sorted_contigs = {cid:cc["contig_length"] for cid, cc in sorted(unassigned_contigs.items(), key=lambda x: x[1]["contig_length"], reverse=True)}
    list_of_contigs_id = list(sorted_contigs.keys())
    list_of_contigs_length = list(sorted_contigs.values())

    # Initialize two pointers
    istart = 0
    jstart = len(list_of_contigs_length)-1
    prev_jstart = jstart + 1
    subset_of_contigs = defaultdict(dict)
    list_of_partitoned_contigs = [] # for valication purpose

    while istart <= jstart:
        # Add one long contig to current chunk
        curr_chunk_length = list_of_contigs_length[istart]
        # Iteratively add smaller contigs until exceed chunk size
        while curr_chunk_length + sum(list_of_contigs_length[jstart:prev_jstart]) <= chunk_size and istart < jstart:
            jstart = jstart - 1

        # Collect all the added shorted contigs
        added_clens = list_of_contigs_length[jstart+1:prev_jstart]
        curr_clens = [curr_chunk_length] + added_clens
        curr_chunk_length += sum(added_clens)

        # Record the list of contig_ids assigned to current chunk_id
        curr_cids = [list_of_contigs_id[istart]] + list_of_contigs_id[jstart+1:prev_jstart]
        subset_of_contigs[chunk_id] = {
            "chunk_id": chunk_id,
            "contigs_id": curr_cids,
            "chunk_length": curr_chunk_length,
            "list_of_contigs_length": curr_clens
        }

        # Update the pointer and chunk_id
        list_of_partitoned_contigs = list_of_partitoned_contigs + curr_clens
        istart = istart + 1
        prev_jstart = jstart + 1
        chunk_id += 1

    assert len(list_of_partitoned_contigs) == len(list_of_contigs_length)
    assert set(list_of_partitoned_contigs) == set(list_of_contigs_length)

    return (subset_of_contigs, chunk_id)

def design_merge_snps_chunks(species_id, contigs_file, chunk_size):

    contigs = scan_fasta(contigs_file)

    # Start with full chunks
    chunk_id = 0
    chunks_of_sites = defaultdict(list)
    unassigned_contigs = defaultdict(dict)
    max_contig_length = 0

    for contig in contigs.values():
        contig_length = contig["length"]
        if contig_length > max_contig_length:
            max_contig_length = contig_length

        contig_id = contig["id"]
        # left closed, right open
        if contig_length < chunk_size:
            unassigned_contigs[contig_id] = {"contig_id": contig_id,
                                             "contig_start": 0,
                                             "contig_end": contig_length,
                                             "contig_length": contig_length}
        else:
            number_of_full_chunks = floor(contig_length/chunk_size)
            for ni, ci in enumerate(range(0, contig_length, chunk_size)):
                if ni == number_of_full_chunks - 1: # last full chunk: carry over
                    chunks_of_sites[chunk_id] = [(species_id, chunk_id, contig_id, ci, contig_length)]
                    chunk_id += 1
                    break
                else:
                    chunks_of_sites[chunk_id] = [(species_id, chunk_id, contig_id, ci, ci+chunk_size)]
                    chunk_id += 1

    if unassigned_contigs:
        # Partition unassigned short contigs into subsets
        dict_of_chunks, chunk_id = partition_contigs_into_chunks(unassigned_contigs, chunk_size, chunk_id)
        # Add the partitioned subsets to chunks
        for chunk_dict in dict_of_chunks.values():
            _chunk_id = chunk_dict["chunk_id"]
            chunks_of_sites[_chunk_id].append((species_id, _chunk_id, -1, chunk_dict["contigs_id"]))
        assert chunk_id == _chunk_id+1

    number_of_chunks = chunk_id
    chunks_of_sites[-1] = (species_id, -1, number_of_chunks, max_contig_length)

    return chunks_of_sites

if __name__ == "__main__":

    CLI=argparse.ArgumentParser()
    CLI.add_argument(
    "--fna",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--genome",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--species",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--run_chunk_size",
    type=int,  
    required=True,
    )
    CLI.add_argument(
    "--merge_chunk_size",
    type=int,  
    required=True,
    )

    args = CLI.parse_args()
    fna_path = args.fna
    genome = args.genome
    species = args.species
    run_chunk_size = args.run_chunk_size
    merge_chunk_size = args.merge_chunk_size

    run_snp_chunks_to_cache = design_run_snps_chunks(species, fna_path, run_chunk_size)
    merge_snp_chunks_to_cache = design_merge_snps_chunks(species, fna_path, merge_chunk_size)

    with open(f'run/chunksize.{run_chunk_size}/{species}/{genome}.json', 'w') as stream:
        json.dump(run_snp_chunks_to_cache, stream)

    with open(f'merge/chunksize.{merge_chunk_size}/{species}/{genome}.json', 'w') as stream:
        json.dump(merge_snp_chunks_to_cache, stream)