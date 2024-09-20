import argparse

def compute_chunks_worker(args):

    violation = "Please do not call compute_chunks_worker directly.  Violation"
    assert args.zzz_worker_mode, f"{violation}:  Missing --zzz_worker_mode arg."

    species_id = args.species
    midas_db = MIDAS_DB(args.midasdb_dir, args.midasdb_name)

    repgenome_for_species = midas_db.uhgg.representatives
    assert species_id in repgenome_for_species, f"{violation}: Species {species_id} is not in the database."
    genome_id = repgenome_for_species[species_id]

    if args.chunk_type == "run_snps":
        contigs_fp = midas_db.fetch_file("representative_genome", species_id, genome_id)
        chunks_to_cache = design_run_snps_chunks(species_id, contigs_fp, args.chunk_size)
    if args.chunk_type == "merge_snps":
        contigs_fp = midas_db.fetch_file("representative_genome", species_id, genome_id)
        chunks_to_cache = design_merge_snps_chunks(species_id, contigs_fp, args.chunk_size)

    dest_filename, _ = get_dest_filename(args.chunk_type, species_id, genome_id)

    local_file = midas_db.get_target_layout(dest_filename, False, species_id, genome_id, args.chunk_size)
    with OutputStream(local_file) as stream:
        json.dump(chunks_to_cache, stream)

    dest_file = midas_db.get_target_layout(dest_filename, True, species_id, genome_id, args.chunk_size)
    if args.upload:
        command(f"aws s3 rm {dest_file}")
        upload(local_file, dest_file)






if __name__ == "__main__":
     
    CLI=argparse.ArgumentParser()
    CLI.add_argument(
    "--species",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--repgenome_fna",
    type=str,  
    required=True,
    )
    CLI.add_argument(
    "--chunk_size",
    type=int,  
    required=True,
    )
    args = CLI.parse_args()

    species_id = args.species
    contigs_fp = args.repgenome_fna

    if args.chunk_type == "run_snps":
        chunks_to_cache = design_run_snps_chunks(species_id, contigs_fp, args.chunk_size)