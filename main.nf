#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.marker_set = "phyeco"
params.db_name = params.db_name
params.db_dir = file(params.db_dir)
params.midas_dir = '/wynton/protected/home/sirota/clairedubin/bin/MIDAS'

params.vsearch_cluster_percents = [99, 95, 90, 85, 80, 75]

// Ensure --db_dir and --midas_dir ends with trailing "/" characters
if (!params.db_dir.endsWith("/")){
    params.db_dir_path = params.db_dir.concat("/")
} else {
    params.db_dir_path = params.db_dir
}
if (!params.midas_dir.endsWith("/")){
    params.midas_dir_path = params.midas_dir.concat("/")
} else {
    params.midas_dir_path = params.midas_dir
}

params.genomes_tsv_path = (params.db_dir_path.concat("genomes.tsv"))
params.marker_model_hmm = params.db_dir_path + "markers_models/" + params.marker_set + "/marker_genes.hmm"
params.bin_path = workflow.launchDir + '/bin'

include { ClusterCentroids as ClusterCentroids } from './modules/cluster_centroids' params(
    db_dir_path: params.db_dir_path
)
include { ClusterCentroids as ReClusterCentroids } from './modules/cluster_centroids' params(
    db_dir_path: params.db_dir_path
)

// TODO:
// check that genome IDs are unique
// check that genome IDs are posix compliant
// check that theres one rep genome per species
// add size for grouptuple for efficiency
// make has_ambiguous_bases a common function


workflow {

    genomes = Channel
        .fromPath(params.genomes_tsv_path)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.genome, row.species, row.representative, row.genome_is_representative) } 

    // get just representative genomes
    genomes.filter { r -> (r[3] == "1") }
        .map{ r -> tuple(r[0], r[1]) }
        .set{ rep_genomes }

    AnnotateGenomes(genomes)
    
    GenerateGeneFeatures(
        AnnotateGenomes.out.genome, 
        AnnotateGenomes.out.species, 
        AnnotateGenomes.out.gff)
    
    HMMMarkerSearch(
        AnnotateGenomes.out.genome, 
        AnnotateGenomes.out.species, 
        AnnotateGenomes.out.faa, 
        AnnotateGenomes.out.ffn)

    // TODO: add a comment about what is in HMMMarkerSearch.out tuple
    ParseHMMMarkers(HMMMarkerSearch.out)

    rep_genomes
        .join(ParseHMMMarkers.out)
        .set{ rep_inferred_markers }

    BuildMarkerDB(
        rep_inferred_markers.collect{ r -> r[2]},
        rep_inferred_markers.collect{ r -> r[3]}
    )

    CleanGenes(
        AnnotateGenomes.out.genome,
        AnnotateGenomes.out.species,
        AnnotateGenomes.out.ffn
        )

    CleanGenes.out.groupTuple(by: 1).set {cleaned_genes_by_species}

    CombineCleanedGenes(
        cleaned_genes_by_species
    )

    Channel
        .fromList(params.vsearch_cluster_percents)
        .max()
        .combine(CombineCleanedGenes.out)
        .set{max_cluster_ch}

    // cluster based on highest clustering value
    ClusterCentroids(
       max_cluster_ch
    ).set {max_cluster_output}


    CleanCentroids(
        max_cluster_output.species,
        max_cluster_output.cluster_pct,
        max_cluster_output.centroids_ffn
    )

    // Recluster at lower thresholds
    max_cluster_val = params.vsearch_cluster_percents.max()
    remaining_clusters_list = params.vsearch_cluster_percents.findAll {it != max_cluster_val}
        
    Channel
        .fromList(remaining_clusters_list)
        .combine(max_cluster_output.recluster_input)
        .set{recluster_ch}

    ReClusterCentroids(
       recluster_ch).set {recluster_output}

    recluster_output.parse_centroid_input
        .concat(max_cluster_output.parse_centroid_input)
        .groupTuple(by: 0)
        .set {uclust_to_parse}

    ParseCentroidInfo(uclust_to_parse)

    ParseCentroidInfo.out
    .join(CleanCentroids.out)
    .join(CombineCleanedGenes.out)
    .set{pipeline_input}

    RunPipelineScript(pipeline_input)

}

process AnnotateGenomes {
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

process GenerateGeneFeatures {
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

process HMMMarkerSearch {
    label 'mem_medium'
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


process ParseHMMMarkers {
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

    script:

    """
    #! /usr/bin/env bash
    set -e
    set -x

    python3 ${params.bin_path}/infer_markers.py --genome ${genome} --species ${species} --hmmsearch_file ${hmmsearch} --annotation_ffn ${ffn}
    
    """

}

process BuildMarkerDB {
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
    #! /usr/bin/env bash
    set -e
    set -x

    cat *.markers.fa >> phyeco.fa
    cat *.markers.fa >> phyeco.map

    hs-blastn index phyeco.fa
    """
}


process CleanGenes {
    label 'mem_low_single_cpu'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/mtest'

    input:
    val(genome)
    val(species)
    path(ffn)

    output:
    tuple val(genome), val(species), path("${genome}.genes.ffn"), path("${genome}.genes.len")

"""
#!/usr/bin/env python3

import Bio.SeqIO

def has_ambiguous_bases(sequence):
    # Check if sequence contains lower-case letters, which usually indicate soft-masked bases
    ambiguous_bases = ['N', 'X', 'n', 'x']
    return any(base in ambiguous_bases for base in sequence)


genome_id = "${genome}"
output_genes = "${genome}.genes.ffn"
output_len = "${genome}.genes.len"
with open(output_genes, 'w') as o_genes, \
        open(output_len, 'w') as o_info:        
    
    for rec in Bio.SeqIO.parse("${ffn}", 'fasta'):
        gene_id = rec.id
        gene_seq = str(rec.seq).upper()
        gene_len = len(gene_seq)
        if gene_len <= 200 or has_ambiguous_bases(gene_seq) or gene_id == '' or gene_id == '|':
            pass
        else:
            o_genes.write(f">{gene_id}\\n{gene_seq}\\n")
            o_info.write(f"{gene_id}\\t{genome_id}\\t{gene_len}\\n")

"""
}

process CombineCleanedGenes {
    label 'mem_low_single_cpu'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/mtest'
    publishDir "${params.db_dir_path}/pangenomes/${species}/temp/vsearch"

    input:
    tuple val(genome), val(species), path('*.genes.ffn'), path('*.genes.len')

    output:
    tuple val(species), path("genes.ffn"), path("genes.len")
    // val(species), emit: species
    // path("genes.ffn"), emit: genes_ffn
    // path("genes.len"), emit: genes_len


    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    cat *.genes.ffn >> genes.ffn
    cat *.genes.len >> genes.len
    
    """

}

process CleanCentroids {
    label 'mem_low_single_cpu'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/mtest'
    publishDir "${params.db_dir_path}/pangenomes/${species}/temp/vsearch/"

    input:
    val(species)
    val(cluster_pct)
    path(centroids_ffn)

    output:
    tuple val(species), path("centroids.${cluster_pct}.clean.ffn"), path("centroids.${cluster_pct}.ambiguous.ffn")

"""
#!/usr/bin/env python3

import Bio.SeqIO

def has_ambiguous_bases(sequence):
    # Check if sequence contains lower-case letters, which usually indicate soft-masked bases
    ambiguous_bases = ['N', 'X', 'n', 'x']
    return any(base in ambiguous_bases for base in sequence)

ffn_in = "${centroids_ffn}"

output_ambiguous = ffn_in.replace('.ffn', '.ambiguous.ffn')
output_clean = ffn_in.replace('.ffn', '.clean.ffn')
with open(output_ambiguous, 'w') as o_ambiguous, \
        open(output_clean, 'w') as o_clean:
    for rec in Bio.SeqIO.parse(ffn_in, 'fasta'):
        c_id = rec.id
        c_seq = str(rec.seq).upper()
        if has_ambiguous_bases(c_seq):
            o_ambiguous.write(f">{c_id}\\n{c_seq}\\n")
        else:
            o_clean.write(f">{c_id}\\n{c_seq}\\n")

"""
}


process ParseCentroidInfo {

    label 'mem_low_single_cpu'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/MIDASv3'
    publishDir "${params.db_dir_path}/pangenomes/${species}/temp/vsearch/"

    input:
    tuple val(species), path('uclust.*.txt')

    output:
    tuple val(species), path("gene_info.txt")

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    python3 ${params.bin_path}/parse_centroids.py uclust.*.txt
    """
}



process RunPipelineScript {

    label 'mem_medium'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/mtest'
    publishDir "${params.db_dir_path}/pangenomes/${species}/temp/cdhit/"
    
    input:
    tuple val(species), path(gene_info), path(clean_centroids), path(amb_centroids), path(genes_ffn), path(genes_len)

    // output:
    // // all of these as a tuple
    // val(genome), emit: genome
    // val(species), emit: species
    // path("${genome}.hmmsearch"), emit: hmmsearch
    // path(ffn)

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x
    
    bash ${params.bin_path}/pipeline.sh \
        ${species} \
        ${gene_info} \
        ${clean_centroids} \
        ${amb_centroids} \
        ${genes_ffn} \
        ${genes_len} \
        ${task.cpus} \
        ${task.memory.toMega()} \
        "${params.midas_dir_path}bin" \

    """



}