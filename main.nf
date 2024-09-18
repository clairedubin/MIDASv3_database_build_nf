#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.marker_set = "phyeco"
params.db_name = params.db_name
params.db_dir = file(params.db_dir)
params.midas_dir = '/wynton/protected/home/sirota/clairedubin/bin/MIDAS'
params.eggnog_db_dir = '/wynton/group/sirota/clairedubin/eggnog'
params.genomad_db_dir="/wynton/protected/home/sirota/clairedubin/databases/genomad_db"
params.resfinder_db_dir="/wynton/protected/home/sirota/clairedubin/databases/resfinder_dbs"
params.resfinder_env_path="/wynton/protected/home/sirota/clairedubin/envs/resfinder_env"

params.vsearch_cluster_percents = [99, 95, 90, 85, 80, 75]
params.max_cluster_val = params.vsearch_cluster_percents.max()


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
include { ClusterCentroids as ClusterCentroidsLowerThresholds } from './modules/cluster_centroids' params(
    db_dir_path: params.db_dir_path
)
include { ReClusterCentroids as ReClusterCentroids } from './modules/cluster_centroids' params(
    db_dir_path: params.db_dir_path
)
include { ReClusterCentroids as ReClusterCentroidsLowerThresholds } from './modules/cluster_centroids' params(
    db_dir_path: params.db_dir_path
)


// TODO:
// check that genome IDs are unique
// check that genome IDs are posix compliant
// check that theres one rep genome per species
// add size for grouptuple for efficiency
// make has_ambiguous_bases a common function
// add outputs for pipeline.sh
// edit directories/outputs in pipeline.sh for speed/redundancy
// check that output files are not empty
//fix hard coded git paths in RunResFinder




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
        .join(ParseHMMMarkers.out, by: 1)
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


    CleanGenes.out.groupTuple(by: 1).set{cleaned_genes_by_species}

    CombineCleanedGenes(cleaned_genes_by_species)

    Channel.of(params.max_cluster_val)
        .combine(CombineCleanedGenes.out.genes_ffn)
        .set{max_cluster_ch}


    // cluster based on highest clustering value
    ClusterCentroids(max_cluster_ch)
        .set {max_cluster_output}

    CleanCentroids(
        max_cluster_output.centroid_ffn,
        max_cluster_output.cluster_pct
    )

    // Recluster at lower thresholds
    remaining_clusters_list = params.vsearch_cluster_percents.findAll {it != params.max_cluster_val}
        
    Channel
        .fromList(remaining_clusters_list)
        .combine(max_cluster_output.centroid_ffn)
        .set{lower_cluster_ch}

    ClusterCentroidsLowerThresholds(lower_cluster_ch)
        .set {lower_cluster_output}

    lower_cluster_output.uclust
        .concat(max_cluster_output.uclust)
        .groupTuple(by: 0)
        .set {uclust_to_parse}

    ParseCentroidInfo(uclust_to_parse)

    ParseCentroidInfo.out
        .join(CleanCentroids.out)
        .join(CombineCleanedGenes.out.genes_ffn)
        .join(CombineCleanedGenes.out.genes_len)
        .set{pipeline_input}

    RefineClusters(pipeline_input)

    // Recluster at lower thresholds        
    Channel
        .fromList(remaining_clusters_list)
        .combine(RefineClusters.out.centroid_ffn)
        .set{lower_cluster_refine_ch}

    ReClusterCentroidsLowerThresholds(lower_cluster_refine_ch)

    ParseHMMMarkers.out
        .map{ r -> tuple(r[0], r[3])}
        .groupTuple(by: 0)
        .set{markers_per_species}

    RefineClusters.out.gene_files
        .combine(ReClusterCentroidsLowerThresholds.out.uclust.groupTuple(by: 0), by: 0)
        .combine(markers_per_species, by: 0)
        .set{recluster_to_parse}

    ParseReclusteredCentroidInfo(recluster_to_parse)

    AugmentPangenomes(ParseReclusteredCentroidInfo.out)
        
    CalculateContigLength(
        AnnotateGenomes.out.fna_tuple.groupTuple(by: 1).map{ r -> tuple(r[1], r[2])})

    RunEggNog(RefineClusters.out.max_centroid_ffn)
    
    RunGeNomad(AnnotateGenomes.out.fna_tuple)
    RunMEFinder(AnnotateGenomes.out.fna_tuple)
    RunResFinder(AnnotateGenomes.out.fna_tuple)

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
    tuple val(genome), val(species), path("${genome}.fna"), emit: fna_tuple
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
    tuple val(genome), val(species), path("${genome}.genes"), emit: genes
    path("${gff}.db")
    

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
    publishDir "${params.db_dir_path}/markers/${params.marker_set}/temp/${species}/${genome}", mode: "copy"
    
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
    publishDir "${params.db_dir_path}/markers/${params.marker_set}/temp/${species}/${genome}", mode: "copy"

    input:
    //make tuple
    val(genome)
    val(species)
    path(hmmsearch)
    path(ffn)

    output:
    tuple val(species), val(genome), path("${genome}.markers.fa"), path("${genome}.markers.map") 

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
    publishDir "${params.db_dir_path}/markers/${params.marker_set}", mode: "copy"

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

    cat *.markers.fa > phyeco.fa
    cat *.markers.fa > phyeco.map

    hs-blastn index phyeco.fa
    """
}


process CleanGenes {
    label 'mem_low_single_cpu'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/mtest'
    publishDir "${params.db_dir_path}/gene_annotations/${species}/${genome}", mode: "copy"


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
    publishDir "${params.db_dir_path}/pangenomes/${species}/temp/vsearch", mode: "copy"

    input:
    tuple val(genome), val(species), path(gene_ffns), path(gene_lens)

    output:
    tuple val(species), path("genes.ffn"), emit: genes_ffn
    tuple val(species), path("genes.len"), emit: genes_len


    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    cat ${gene_ffns} > "genes.ffn"
    cat ${gene_lens} > "genes.len"
    
    if [ ! -s 'genes.ffn' ]; then
        echo ERROR: genes.ffn is empty
        exit 1 
    elif [ ! -s 'genes.len' ]; then
        echo ERROR: genes.len is empty
        exit 1 

    fi
    """

}

process CleanCentroids {
    label 'mem_low_single_cpu'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/mtest'
    publishDir "${params.db_dir_path}/pangenomes/${species}/temp/vsearch/", mode: "copy"

    input:
    tuple val(species), path(centroid_ffn)
    val(cluster_pct)


    output:
    tuple val(species), path("centroids.${cluster_pct}.clean.ffn"), path("centroids.${cluster_pct}.ambiguous.ffn")

"""
#!/usr/bin/env python3

import Bio.SeqIO
import os

def has_ambiguous_bases(sequence):
    # Check if sequence contains lower-case letters, which usually indicate soft-masked bases
    ambiguous_bases = ['N', 'X', 'n', 'x']
    return any(base in ambiguous_bases for base in sequence)

ffn_in = "${centroid_ffn}"

assert os.path.getsize(ffn_in) > 0

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
    publishDir "${params.db_dir_path}/pangenomes/${species}/temp/vsearch/", mode: "copy"

    input:
    tuple val(species), path(uclust_files)

    output:
    tuple val(species), path("gene_info.txt")

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    python3 ${params.bin_path}/parse_centroids.py ${uclust_files}

    """
}



process RefineClusters {

    label 'mem_medium'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/mtest'
    publishDir (
        path: "${params.db_dir_path}/pangenomes/${species}/",
        mode: "copy",
        saveAs: { fn ->
            if (fn == "centroids.ffn") { "${fn}" }
            else { "temp/cdhit/${fn}" }
        }
    )
    
    input:
    tuple val(species), path(gene_info), path(clean_centroids), path(amb_centroids), path(genes_ffn), path(genes_len)

    output:
    tuple val(species), path("centroids.*.refined.ffn"), emit: centroid_ffn 
    tuple val(species), path("genes.len.txt"), path("gene_info.txt"), emit: gene_files
    tuple val(species), path("centroids.ffn"), emit: max_centroid_ffn 
    path("genes.ffn")

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


    if [ ! -e PIPELINE_SUCCESS ]; then
        echo ERROR: Cluster refinement not complete
        exit 1 
    fi

    cp "centroids.${params.max_cluster_val}.refined.ffn" "centroids.ffn"

    """
}

process ParseReclusteredCentroidInfo {

    label 'mem_low_single_cpu'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/MIDASv3'
    publishDir "${params.db_dir_path}/pangenomes/${species}/", mode: "copy"

    input:
    tuple val(species), path(genes_len), path(genes_info), path(uclust_files), path(marker_map_files)

    output:
    tuple val(species), path('genes_info.tsv')

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    python3 ${params.bin_path}/parse_reclustered_centroids.py \
    --gene_info_file ${genes_info} \
    --max_percent ${params.max_cluster_val} \
    --gene_length_file ${genes_len} \
    --uclust_files ${uclust_files} \
    --genome_marker_files ${marker_map_files}

    """
}

process AugmentPangenomes {

    label 'mem_low_single_cpu'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/MIDASv3'
    publishDir "${params.db_dir_path}/pangenomes/${species}/augment/", mode: "copy"

    input:
    tuple val(species), path(gene_info_tsv), path(fna_list)

    output:
    tuple val(species), path('clusters_{xx}_info.tsv')

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    cluster_pct_list=\$(echo ${params.vsearch_cluster_percents} | tr -d '[],')

    python3 ${params.bin_path}/augment_pangenome.py \
    --gene_info_file ${gene_info_tsv} \
    --cluster_thresholds \${cluster_pct_list} \
    """
}

process CalculateContigLength {

    label 'mem_low_single_cpu'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/MIDASv3'
    publishDir "${params.db_dir_path}/pangenomes/${species}/", mode: "copy"

    input:
    tuple val(species), path(fna_list)

    output:
    tuple val(species), path("contigs.len")

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    python3 ${params.bin_path}/calculate_contig_length.py \
    --fna_list ${fna_list} \
    """
}


process RunEggNog {

    label 'mem_medium'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/eggnog'
    publishDir "${params.db_dir_path}/pangenomes_annotation/02_eggnog/${species}", mode: "copy"
    
    input:
    tuple val(species), path(max_centroid_ffn)

    output:
    val(species), path("${species_id}.emapper.annotations")

    script:
    """
    #! /usr/bin/env bash
    set -e
    
    emapper.py \
    -i ${max_centroid_ffn} --itype CDS \
    -m diamond --sensmode more-sensitive \
    --data_dir ${params.eggnog_db_dir} \
    --output ${species} --override \
    --dbmem --pfam_realign realign \
    --cpu ${task.cpus} \

    """
}

process RunGeNomad {
    
    label 'mem_medium'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/genomad'
    publishDir "${params.db_dir_path}/pangenomes_annotations/01_mge/${species}/${genome}/genomad_output", mode: "copy"
    
    input:
    tuple val(genome), val(species), path(fna)

    output:
    tuple val(genome), val(species), \
    path("${genome}_summary/GCA_007135585.1_virus_genes.tsv"), 
    path("${genome}_summary/GCA_007135585.1_plasmid_genes.tsv"), emit: genomad_gene_files

    path("${genome}_*")

    script:
    """
    #! /usr/bin/env bash
    set -e
    
    genomad end-to-end --threads ${task.cpus} --cleanup --enable-score-calibration ${fna} '.' ${params.genomad_db_dir} 

    """
}

process RunMEFinder {
    
    label 'mem_medium'
    errorStrategy 'finish'
    publishDir "${params.db_dir_path}/pangenomes_annotations/01_mge/${species}/${genome}/mefinder_output", mode: "copy"
    
    input:
    tuple val(genome), val(species), path(fna)

    output:
    tuple val(genome), val(species), path("mefinder.csv")

    script:
    """
    #! /usr/bin/env bash
    set -e
    module load CBI blast
    source ${params.resfinder_env_path}/bin/activate
    
    mefinder find --contig ${fna} -t ${task.cpus} mefinder

    """
}

process RunResFinder {
    
    label 'mem_medium'
    errorStrategy 'finish'
    publishDir "${params.db_dir_path}/pangenomes_annotations/01_mge/${species}/${genome}/resfinder_output", mode: "copy"
    
    input:
    tuple val(genome), val(species), path(fna)

    output:
    tuple val(genome), val(species), path('ResFinder_results_tab.txt')

    script:
    """
    #! /usr/bin/env bash
    set -e
    export GIT_PYTHON_GIT_EXECUTABLE="/wynton/protected/home/sirota/clairedubin/bin/git-2.39.5/git"
    export PATH=$PATH:/wynton/protected/home/sirota/clairedubin/bin/git-2.39.5
    source ${params.resfinder_env_path}/bin/activate
    
    python -m resfinder -ifa ${fna} -o resfinder_output \
    -s Other -l 0.6 -t 0.8 --acquired -db_res ${params.resfinder_db_dir}/resfinder_db \
    -d -db_disinf ${params.resfinder_db_dir}/disinfinder_db

    """
}


process ParsePangenomeAnnotations {

    label 'mem_low_single_cpu'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/MIDASv3'
    // publishDir "${params.db_dir_path}/pangenomes_annotation/03_processed/${species}/", mode: "copy"
    //     "panannot_tempfile":             f"pangenomes_annotation/03_processed/{species_id}/{component}/{genome_id}",
    // // publishDir (
    // //     path: "${params.db_dir_path}/pangenomes_annotation/03_processed/${species}/",
    // //     mode: "copy",
    // //     saveAs: { fn ->
    // //         if (fn.startsWith("A.tab")) { "genomad_virus/${genome}/${fn}" }
    // //         if (fn.startsWith("A.tab")) { "genomad_virus/${genome}/${fn}" }
    // //         if (fn.startsWith("A.tab")) { "genomad_virus/${genome}/${fn}" }
    // //         if (fn.startsWith("A.tab")) { "genomad_virus/${genome}/${fn}" }
    // //         if (fn.startsWith("A.tab")) { "genomad_virus/${genome}/${fn}" }
    // //         else { "unassigned/${fn}" }
    // //     }
    // // )
//        "resfinder_results":             f"pangenomes_annotation/01_mge/{species_id}/{genome_id}/resfinder_output/ResFinder_results_tab.txt",
        // "mefinder_results":              f"pangenomes_annotation/01_mge/{species_id}/{genome_id}/mefinder_output/mefinder.csv",
//         "genomad_virus_genes":           f"pangenomes_annotation/01_mge/{species_id}/{genome_id}/genomad_output/{genome_id}_summary/{genome_id}_virus_genes.tsv",
//         "genomad_plasmid_genes":         f"pangenomes_annotation/01_mge/{species_id}/{genome_id}/genomad_output/{genome_id}_summary/{genome_id}_plasmid_genes.tsv",

    // cleangenes.out: val(genome), val(species), path("${genome}.genes.ffn")
// mefinder.out: genome, species, mefinder
// resfinder.out: genome, species, resfinder
// rungenomad.out.genomad_gene_files: genome, species, virus genes, plasmid genes
    // CalculateContigLength.out: species, "contigs.len"

CleanGenes.out.genes
    .join(RunMEFinder.out, by: [0, 1])
    .join(RunResFinder.out, by: [0, 1])
    .join(RunGeNomad.out, by: [0, 1])
    .map{it.swap(1,0)} //swap order of genome and species
    .join(CalculateContigLength.out, by: 0)
    .join(RunEggNog.out, by: 0)
    .join(AugmentPangenomes.out, by: 0)


    input:
    tuple val(genome), val(species), path(genes_file), \
    path(mefinder_csv), path(resfinder_txt), \
    path(genomad_virus), path(genomad_plasmid), \
    path(contig_len_file), path(eggnog_results), path(max_cluster_info)

    output:
    tuple val(species), 

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    python3 ${params.bin_path}/parse_pangenome_annotations.py \
    --species ${species} \
    --genome ${genome} \
    --genes_file ${genes_file} \
    --contig_len_file ${contig_len_file} \
    --mefinder_csv ${mefinder_csv} \
    --resfinder_txt ${resfinder_txt} \
    --genomad_virus_file ${genomad_virus} \
    --genomad_plasmid_file ${genomad_plasmid} \
    --max_cluster_info_file ${max_cluster_info}
    --eggnog_results_file ${eggnog_results}
    """
}



    

    // --eggnog_results_file ${eggnog_results} 
    // --max_cluster_info_file ${max_cluster_info} \


// process ParseEggNogAnnotations {
//     label 'mem_low_single_cpu'
//     errorStrategy 'finish'
//     conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/mtest'
//     publishDir "${params.db_dir_path}/pangenomes/${species}/temp/vsearch/", mode: "copy"

//     input:
//     tuple val(genome), val(species), path(max_cluster_info_file), \
//     path(contig_len_file), path(genes_file), path(eggnog_results_file)

//     output:
//     tuple val(species), path("eggnog/${genome}")

// """
// #!/usr/bin/env python3

// import pandas as pd

// genome_id = ${genome}
// species_id = ${species}
// cluster_info_fp = ${max_cluster_info_file}
// contig_len_fp = ${contig_len_file}
// gene_feature_fp = ${genes_file}
// eggnog_results_file = ${eggnog_results_file}

// contig_len = pd.read_csv(contig_len_fp, sep='\t')
// gene_features = pd.read_csv(gene_feature_fp, sep='\t')

// centroids_99 = pd.read_csv(cluster_info_fp, sep='\t', usecols=[0]) # Only keep the list of centroid_99
// centroids_99 = pd.merge(centroids_99, gene_features, left_on='centroid_99', right_on='gene_id', how='inner')
// centroids_99 = pd.merge(centroids_99, contig_len[['contig_id', 'contig_length']], left_on='contig_id', right_on='contig_id', how='inner')
// centroids_99 = centroids_99[['contig_id', 'start', 'end', 'centroid_99', 'strand', 'gene_type', 'contig_length']]

// df = scan_eggnog(eggnog_results_file)
// local_dest = f"eggnog/{genome_id}"
// merged_df  = df.merge(centroids_99, left_on='#query', right_on='centroid_99', how='inner')
// merged_df = merged_df.drop(columns=['centroid_99'])
// merged_df.to_csv(local_dest, sep='\t', index=False, header=False)

// """
// }