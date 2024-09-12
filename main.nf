#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.marker_set = "phyeco"
params.db_name = params.db_name
params.db_dir = file(params.db_dir)


// Ensure --db_dir ends with trailing "/" characters
if (!params.db_dir.endsWith("/")){
    params.db_dir_path = params.db_dir.concat("/")
} else {
    params.db_dir_path = params.db_dir
}

params.genomes_tsv_path = (params.db_dir_path.concat("genomes.tsv"))
params.marker_model_hmm = params.db_dir_path + "markers_models/" + params.marker_set + "/marker_genes.hmm"
params.bin_path = workflow.launchDir + '/bin'

// //TODO:
// check that genome IDs are unique
// check that genome IDs are posix compliant
// check that theres one rep genome per species


workflow {

    genomes = Channel
        .fromPath(params.genomes_tsv_path)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.genome, row.species, row.representative, row.genome_is_representative) } 

    // get just representative genomes
    genomes.filter { r -> (r[3] == "1") }
        .map{ r -> tuple(r[0], r[1]) }
        .set{ rep_genomes }

    annotateGenomes(genomes)
    
    generateGeneFeatures(
        annotateGenomes.out.genome, 
        annotateGenomes.out.species, 
        annotateGenomes.out.gff)
    
    hmmMarkerSearch(
        annotateGenomes.out.genome, 
        annotateGenomes.out.species, 
        annotateGenomes.out.faa, 
        annotateGenomes.out.ffn)

    // TODO: add a comment about what is in hmmMarkerSearch.out tuple
    inferMarkers(hmmMarkerSearch.out)

    rep_genomes
        .join(inferMarkers.out)
        .set{ rep_inferred_markers }

    rep_inferred_markers.view()

    buildMarkerDB(
        rep_inferred_markers.collect{ r -> r[2]},
        rep_inferred_markers.collect{ r -> r[3]}
    )

}

process annotateGenomes {
    label 'mem_low'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/mtest'
    publishDir "${params.db_dir_path}/gene_annotations/${species}/${genome}", mode: "copy"

    input:
    tuple val(genome), val(species), val(representative), val(genome_is_representative)

    output:
    /// can also make a big tuple with everything i need later on
    // tuple val(genome), val(species), emit: input_annotation
    val(genome), emit: genome
    val(species), emit: species
    path("${genome}.gff"), emit: gff
    path("${genome}.faa"), emit: faa
    path("${genome}.ffn"), emit: ffn
    path("${genome}.tsv"), emit: tsv
    path("${genome}.fna"), emit: fna
    path("${genome}.log")

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    prokka --kingdom Bacteria --metagenome --cpus ${task.cpus} --prefix ${genome} --locustag ${genome} --outdir . --compliant --force "${params.db_dir_path}/cleaned_imports/${species}/${genome}/${genome}.fasta"

    """
}

process generateGeneFeatures {
    label 'mem_low_single_cpu'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/mtest'
    publishDir "${params.db_dir_path}/gene_annotations/${species}/${genome}", mode: "copy"

    input:
    val(genome)
    val(species)
    path(gff)

    output:
    //make tuple: genome, species, gffdb, genomegenes
    path("${gff}.db")
    path("${genome}.genes")

"""
#!/usr/bin/env python3

import gffutils

db = gffutils.create_db("${gff}", "${gff}.db")
to_write = "\\t".join(["gene_id", "contig_id", "start", "end", "strand", "gene_type"]) + "\\n"
for feature in db.all_features():
    if feature.source == "prokka":
        continue
    if "ID" not in feature.attributes: #CRISPR
        continue
    seqid = feature.seqid
    start = feature.start
    stop = feature.stop
    strand = feature.strand
    gene_id = feature.attributes['ID'][0]
    locus_tag = feature.attributes['locus_tag'][0]
    assert gene_id == locus_tag
    gene_type = feature.featuretype
    to_write += "\\t".join([gene_id, seqid, str(start), str(stop), strand, gene_type]) + "\\n"

with open("${genome}.genes", 'w') as f:
    f.write(to_write)

"""
 }

process hmmMarkerSearch {
    label 'mem_med'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/mtest'
    publishDir "${params.db_dir_path}/markers/${params.marker_set}/temp/${species}/${genome}"
    
    input:
    val(genome)
    val(species)
    path(faa)
    path(ffn)

    output:
    // all of these as a tuple
    val(genome), emit: genome
    val(species), emit: species
    path("${genome}.hmmsearch"), emit: hmmsearch
    path(ffn)

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x
    
    hmmsearch --noali --cpu ${task.cpus} --domtblout "${genome}.hmmsearch" ${params.marker_model_hmm} ${faa}

    """
}


process inferMarkers {
    label 'mem_low_single_cpu'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/MIDASv3'
    publishDir "${params.db_dir_path}/markers/${params.marker_set}/temp/${species}/${genome}"

    input:
    //make tuple
    val(genome)
    val(species)
    path(hmmsearch)
    path(ffn)

    output:
    tuple val(genome), path("${genome}.markers.fa"), path("${genome}.markers.map")
    // val(genome), emit: genome
    // val(species), emit: species
    // path("${genome}.markers.fa"), emit: markers_fa
    // path("${genome}.markers.map"), emit: markers_map

    script:

    """
    #! /usr/bin/env bash
    set -e
    set -x

    python3 ${params.bin_path}/infer_markers.py --genome ${genome} --species ${species} --hmmsearch_file ${hmmsearch} --annotation_ffn ${ffn}
    
    """

}

process buildMarkerDB {
    label 'mem_low_single_cpu'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/mtest'
    publishDir "${params.db_dir_path}/markers/${params.marker_set}"

    input:
    path '*.markers.fa'
    path '*.markers.map'

    output:
    path("phyeco.fa")
    path("phyeco.map")


    script:
    """
    cat *.markers.fa >> phyeco.fa
    cat *.markers.fa >> phyeco.map

    hs-blastn index phyeco.fa
    """
}
