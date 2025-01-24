#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// TODO:
// add size for grouptuple for efficiency
// containerize stuff
// clean up files in bin 
// change eggnog db loading for efficiency
// make prokka accept either .fa or .fasta

//// PARAMS
params.centroid_cluster_percents = [99, 95, 90, 85, 80, 75]
params.max_cluster_val = params.centroid_cluster_percents.max()
params.run_chunk_size = 1000000
params.merge_chunk_size = 500000
params.resfinder_min_cov = 0.6
params.resfinder_identity_threshold = 0.8
params.hmmsearch_min_cov = 0.00
params.hmmsearch_max_evalue = 1e-5

//// PATH PARAMS
params.marker_set = "phyeco"
params.db_name = "localdb"
params.db_dir = "/wynton/protected/home/sirota/clairedubin/MIDAS3_nextflow/test_db"
params.midas_dir = "/wynton/protected/home/sirota/clairedubin/bin/MIDAS"

// conda env for most processes
params.conda_env_path = "/wynton/protected/home/sirota/clairedubin/anaconda3/envs/MIDASv3"

// eggnog has its own conda env
params.eggnog_db_dir = "/wynton/group/sirota/clairedubin/eggnog"
params.eggnog_conda_dir = '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/eggnog'

// genomad has its own conda env
params.genomad_db_dir = "/wynton/protected/home/sirota/clairedubin/databases/genomad_db_v1.5/genomad_db"
params.genomad_conda_dir =  '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/genomad'

// resfinder has its own venv
params.resfinder_env_dir = "/wynton/protected/home/sirota/clairedubin/envs/resfinder_env"
params.resfinder_db_dir = "/wynton/protected/home/sirota/clairedubin/databases/resfinder_dbs"

// other necessary programs
params.blastn_dir = "/wynton/protected/home/sirota/clairedubin/bin/ncbi-blast-2.14.1+/bin"
params.git_dir = "/wynton/protected/home/sirota/clairedubin/bin/git-2.39.5"

params.bin_dir = workflow.launchDir + "/bin"

// Ensure directories end with trailing "/" characters
params.findAll { key, _ -> key.endsWith("_dir")}
    .each { key, path ->
    if (!path.endsWith("/")) {
        params[key.replace("_dir", "_path")] = "${path}/"
    } else {
        params[key.replace("_dir", "_path")] = "${path}"
    }
}

// Subfiles within above directories
params.marker_path = params.db_path + "markers_models/" + params.marker_set + '/'
params.marker_model_hmm_path = params.marker_path + "marker_genes.hmm"
params.genomes_tsv_path = params.db_path + "genomes.tsv"
params.git_executable_path = params.git_path + "git"
params.eggnog_dmnd_db_path = params.eggnog_db_path + "eggnog_prokaryote_only.dmnd"

// Check that all paths exist
params.findAll { key, _ -> key.endsWith("_path")}
    .each { key, path ->
    def file = new File(path.toString())
    if (!file.exists()) {
        throw new IllegalStateException("Error: ${key} does not exist at ${path}")
    } 
}

include { check_input } from './modules/check_input' params(
    db_path: params.db_path
)
include { ClusterCentroids as ClusterCentroids } from './modules/cluster_centroids' params(
    db_path: params.db_path,
    conda_env_path: params.conda_env_path
)
include { ClusterCentroids as ClusterCentroidsLowerThresholds } from './modules/cluster_centroids' params(
    db_path: params.db_path,
    conda_env_path: params.conda_env_path

)

workflow {

    //// CHECK GENOMES.TSV FORMATTING 
    check_input(params.genomes_tsv_path)

    genomes = Channel
        .fromPath(params.genomes_tsv_path)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.genome, row.species, row.representative, row.genome_is_representative) } 

    // Filter to only representative genomes for certain steps of pipeline
    genomes.filter { r -> (r[3] == "1") }
        .map{ r -> tuple(r[0], r[1]) } // tuple of (genome, species)
        .set{ rep_genomes }

    //// ANNOTATION ////

    AnnotateGenomes(genomes)
    
    GenerateGeneFeatures(AnnotateGenomes.out.gff_tuple)

    CleanGenes(AnnotateGenomes.out.ffn_tuple)

    CombineCleanedGenes(CleanGenes.out.groupTuple(by: 1))

    RunGeNomad(AnnotateGenomes.out.fna_tuple)
    RunMEFinder(AnnotateGenomes.out.fna_tuple)
    RunResFinder(AnnotateGenomes.out.fna_tuple)

    CalculateContigLength(
        AnnotateGenomes.out.fna_tuple
            .groupTuple(by: 1)
            .map{ r -> tuple(r[1], r[2])}
            )

    //// BUILD MARKER DB ////

    HMMMarkerSearch(AnnotateGenomes.out.faa_ffn_tuple)

    ParseHMMMarkers(HMMMarkerSearch.out)

    rep_genomes
        .join(ParseHMMMarkers.out, by: [0,1])
        .set{ rep_inferred_markers }

    BuildMarkerDB(
        rep_inferred_markers.collect{ r -> r[2]},
        rep_inferred_markers.collect{ r -> r[3]}
    )

    //// CLUSTERING ////

    // Initial clustering based on highest clustering value (usually 99)
    Channel.of(params.max_cluster_val)
        .combine(CombineCleanedGenes.out.genes_ffn)
        .set{max_cluster_ch}
    ClusterCentroids(max_cluster_ch)
        .set {max_cluster_output}

    CleanCentroids(
        max_cluster_output.centroid_ffn,
        max_cluster_output.cluster_pct
    )

    // Clustering C99 output at lower thresholds
    remaining_clusters_list = params.centroid_cluster_percents.findAll {it != params.max_cluster_val}
        
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
        .join(CleanCentroids.out.cleaned_centroids)
        .join(CombineCleanedGenes.out.genes_ffn)
        .join(CombineCleanedGenes.out.genes_len)
        .set{refine_clusters_input}

    RefineClusters(refine_clusters_input)

    // Redoing clustering at lower thresholds        
    Channel
        .fromList(remaining_clusters_list)
        .combine(RefineClusters.out.centroid_ffn)
        .set{lower_cluster_refine_ch}

    ReClusterCentroids(lower_cluster_refine_ch)

    ParseHMMMarkers.out
        .map{ r -> tuple(r[1], r[3])}
        .groupTuple(by: 0)
        .set{markers_per_species}

    RefineClusters.out.gene_files
        .combine(ReClusterCentroids.out.uclust.groupTuple(by: 0), by: 0)
        .combine(markers_per_species, by: 0)
        .set{recluster_to_parse}

    ParseReclusteredCentroidInfo(recluster_to_parse)

    AugmentPangenomes(    
        rep_genomes.map{r -> tuple(r[1], r[0])}
        .combine(ParseReclusteredCentroidInfo.out.reclustered_centroid_info, by: 0)
        )
        
    ////// EGGNOG ANNOTATION OF CENTROIDS ////

    RunEggNog(RefineClusters.out.centroid_ffn)
    
    //// COMBINE ANNOTATION RESULTS ////

    GenerateGeneFeatures.out.genes
        .combine(RunMEFinder.out.mefinder_output, by: [0,1])
        .combine(RunResFinder.out.resfinder_output, by: [0,1])
        .combine(RunGeNomad.out.genomad_gene_files, by: [0,1])
        .map{it.swap(1,0)} // swap genome and species columns to match order of the following inputs
        .combine(CalculateContigLength.out, by: 0)
        .combine(RunEggNog.out.eggnog_output, by: 0)
        .combine(AugmentPangenomes.out.max_cluster_info, by: 0)
        .map{it.swap(1,0)} // swap back genome and species
        .set{annotations_to_parse}

    ParsePangenomeAnnotations(annotations_to_parse)

    CombinePangenomeAnnotations(ParsePangenomeAnnotations.out.groupTuple(by: 0))

    CombinePangenomeAnnotations.out
        .combine(ParseReclusteredCentroidInfo.out.reclustered_centroid_info, by: 0)
        .combine(AugmentPangenomes.out.cluster_info_file_list, by: 0)
        .set{enhance_pangenome_input}
    
    EnhancePangenome(enhance_pangenome_input)

    //// COMPUTE CHUNKS ////

    rep_genomes
        .join(AnnotateGenomes.out.fna_tuple, by: [0,1])
        .set{rep_genome_fna_tuples}

    ComputeRunSNPsChunks(rep_genome_fna_tuples)
    ComputeMergeSNPsChunks(rep_genome_fna_tuples)
    
}

process AnnotateGenomes {

    label 'mem_medium'
    publishDir "${params.db_path}/gene_annotations/${species}/${genome}", mode: "copy"

    input:
    tuple val(genome), val(species), val(representative), val(genome_is_representative)

    output:
    tuple val(genome), val(species), path("${genome}.gff"), emit: gff_tuple
    tuple val(genome), val(species), path("${genome}.ffn"), emit: ffn_tuple
    tuple val(genome), val(species), path("${genome}.fna"), emit: fna_tuple
    tuple val(genome), val(species), path("${genome}.faa"), path("${genome}.ffn"), emit: faa_ffn_tuple
    path("${genome}.tsv")

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    prokka \
    --kingdom Bacteria \
    --metagenome \
    --cpus ${task.cpus} \
    --prefix ${genome} \
    --locustag ${genome} \
    --outdir . \
    --compliant \
    --force \
    "${params.db_path}/cleaned_imports/${species}/${genome}/${genome}.fasta"

    """
}

process GenerateGeneFeatures {
    label 'single_cpu'
    publishDir "${params.db_path}/gene_annotations/${species}/${genome}", mode: "copy"

    input:
    tuple val(genome), val(species), path(gff)

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
    publishDir "${params.db_path}/markers/${params.marker_set}/temp/${species}/${genome}", mode: "copy"
    
    input:
    tuple val(genome), val(species), path(faa), path(ffn)

    output:
    tuple val(genome), val(species), path("${genome}.hmmsearch"), path(ffn)

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x
    
    hmmsearch --noali --cpu ${task.cpus} --domtblout "${genome}.hmmsearch" ${params.marker_model_hmm_path} ${faa}

    """
}

process ParseHMMMarkers {
    label 'single_cpu'
    publishDir "${params.db_path}/markers/${params.marker_set}/temp/${species}/${genome}", mode: "copy"

    input:
    tuple val(genome), val(species), path(hmmsearch), path(ffn)

    output:
    tuple val(genome), val(species), path("${genome}.markers.fa"), path("${genome}.markers.map") 

    script:

    """
    #! /usr/bin/env bash
    set -e
    set -x

    python3 ${params.bin_path}/infer_markers.py \
    --genome ${genome} \
    --species ${species} \
    --hmmsearch_file ${hmmsearch} \
    --annotation_ffn ${ffn} \
    --hmmsearch_max_evalue ${params.hmmsearch_max_evalue} \
    --hmmsearch_min_cov ${params.hmmsearch_min_cov}

    """

}

process BuildMarkerDB {
    label 'single_cpu'
    publishDir "${params.db_path}/markers/${params.marker_set}", mode: "copy"

    input:
    path '*.markers.fa'
    path '*.markers.map'

    output:
    path("phyeco.*")

    script:

    """
    #! /usr/bin/env bash
    set -e
    set -x

    cat *.markers.fa > phyeco.fa
    cat *.markers.map > phyeco.map

    hs-blastn index phyeco.fa
    """
}

process CleanGenes {

    label 'single_cpu'
    publishDir "${params.db_path}/gene_annotations/${species}/${genome}", mode: "copy"

    input:
    tuple val(genome), val(species), path(ffn)

    output:
    tuple val(genome), val(species), path("${genome}.genes.ffn"), path("${genome}.genes.len")


"""
#!/usr/bin/env python3

import Bio.SeqIO
import sys
sys.path.append("${params.bin_path}")
from common import has_ambiguous_bases

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
    label 'single_cpu'
    publishDir "${params.db_path}/pangenomes/${species}/vsearch", mode: "copy"

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

    cat ${gene_ffns} > genes.ffn
    cat ${gene_lens} > genes.len
    
    """
}

process CleanCentroids {
    label 'single_cpu'
    publishDir "${params.db_path}/pangenomes/${species}/temp/vsearch/", mode: "copy"

    input:
    tuple val(species), path(centroid_ffn)
    val(cluster_pct)


    output:
    tuple val(species), \
    path("centroids.${cluster_pct}.clean.ffn"), \
    path("centroids.${cluster_pct}.ambiguous.ffn"), emit: cleaned_centroids

"""
#!/usr/bin/env python3

import Bio.SeqIO
import os
import sys
sys.path.append("${params.bin_path}")
from common import has_ambiguous_bases

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

    label 'single_cpu'
    publishDir "${params.db_path}/pangenomes/${species}/temp/vsearch/", mode: "copy"

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
    publishDir (
        path: "${params.db_path}/pangenomes/${species}",
        mode: "copy",
    )
    
    input:
    tuple val(species), path(gene_info), path(clean_centroids), path(amb_centroids), path(genes_ffn), path(genes_len)

    output:
    tuple val(species), path("centroids.ffn"), emit: centroid_ffn 
    tuple val(species), path("temp/cdhit/genes.len"), path("temp/cdhit/gene_info.txt"), emit: gene_files
    path("temp/cdhit/*")
    path("temp/cdhit/centroids.${params.max_cluster_val}.ffn")

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
        "${params.midas_path}bin" 

    cp "temp/cdhit/centroids.${params.max_cluster_val}.ffn" centroids.ffn


    """
}

process ReClusterCentroids {
    label 'mem_high'
    publishDir "${params.db_path}/pangenomes/${species}/temp/", mode: "copy"

    input:
    tuple val(cluster_pct), val(species), path(genes_ffn)


    output:
    tuple val(species), path("centroids.${cluster_pct}.ffn"), emit: centroid_ffn
    tuple val(species), path("uclust.${cluster_pct}.txt"), emit: uclust


    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    cluster_prop="\$(awk "BEGIN {printf \\"%.2f\\", ${cluster_pct} / 100}")"
    vsearch --cluster_fast ${genes_ffn} --id \${cluster_prop} --threads ${task.cpus} --centroids centroids.${cluster_pct}.ffn --uc uclust.${cluster_pct}.txt

    """
}

process ParseReclusteredCentroidInfo {

    label 'single_cpu'
    publishDir "${params.db_path}/pangenomes/${species}/", mode: "copy"

    input:
    tuple val(species), path(genes_len), path(genes_info), path(uclust_files), path(marker_map_files)

    output:
    tuple val(species), path('genes_info.tsv'), emit: reclustered_centroid_info
    path('temp/*')

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    echo ""
    echo ""
    echo ""
    echo ""


    python3 ${params.bin_path}/parse_reclustered_centroids.py \
    --gene_info_file ${genes_info} \
    --max_percent ${params.max_cluster_val} \
    --gene_length_file ${genes_len} \
    --uclust_files ${uclust_files} \
    --genome_marker_files ${marker_map_files}

    mkdir temp
    cp genes_info.tsv temp/genes_info.tsv
    mv mapfile temp
    cp temp/mapfile temp/markers.map

    """
}

process AugmentPangenomes {

    label 'single_cpu'
    publishDir "${params.db_path}/pangenomes/${species}/augment/", mode: "copy"

    input:
    tuple val(species), val(rep_genome), path(gene_info_tsv)

    output:
    tuple val(species), path("clusters_${params.max_cluster_val}_info.tsv"), emit: max_cluster_info
    tuple val(species), path('clusters_*_info.tsv'), emit: cluster_info_file_list
    path("centroid_*_stacked_matrix.tsv")
    path("centroid_${params.max_cluster_val}_lifted_to_${rep_genome}.tsv")

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    cluster_pct_list=\$(echo ${params.centroid_cluster_percents} | tr -d '[],')

    python3 ${params.bin_path}/augment_pangenome.py \
    --gene_info_file ${gene_info_tsv} \
    --cluster_thresholds \${cluster_pct_list} \
    --rep_genome ${rep_genome}
    """
}

process CalculateContigLength {

    label 'single_cpu'
    publishDir "${params.db_path}/pangenomes/${species}/", mode: "copy"

    input:
    tuple val(species), path(fna_list)

    output:
    tuple val(species), path("contigs.len")

"""
#!/usr/bin/env python3
from Bio import SeqIO

fna_list = "${fna_list}".split(' ')
dict_of_contig_length = {}

for genome_fna in fna_list:
    clen = {}
    genome_id = genome_fna.replace('.fna','')
    for rec in SeqIO.parse(genome_fna, 'fasta'):
        clen[rec.id] = len(rec.seq)
    dict_of_contig_length[genome_id] = clen

with open("contigs.len", 'w') as f:
        f.write("\\t".join(["genome_id", "contig_id", "contig_length"]) + "\\n")
        for gid, r in dict_of_contig_length.items():
            for cid, clen in r.items():
                vals = [gid, cid, str(clen)]
                f.write("\\t".join(vals) + "\\n")
"""
}

process RunEggNog {

    ///this takes a lot of time, mostly loading the annotation db into memory
    ///should explore other ways to just load db once or twice
    //also note i hard coded scratch dir below for wynton compute nodes

    label 'mem_very_high'
    conda "${params.eggnog_conda_path}"
    publishDir "${params.db_path}/pangenomes_annotation/02_eggnog/${species}", mode: "copy"
    
    input:
    tuple val(species), path(centroid_ffn)

    output:
    tuple val(species), path("${species}.emapper.annotations"), emit: eggnog_output
    path("${species}.emapper*")

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    time emapper.py \
    -i ${centroid_ffn} --itype CDS \
    -m diamond --sensmode more-sensitive \
    --data_dir ${params.eggnog_db_path} \
    --output ${species} --override \
    --pfam_realign realign \
    --cpu ${task.cpus} \
    --dmnd_db ${params.eggnog_dmnd_db_path}

    """
}


process RunGeNomad {
    
    label 'mem_high'
    conda "${params.genomad_conda_path}"
    publishDir "${params.db_path}/pangenomes_annotation/01_mge/${species}/${genome}/genomad_output", mode: "copy"
    
    input:
    tuple val(genome), val(species), path(fna)

    output:
    tuple val(genome), val(species), \
    path("${genome}_summary/${genome}_virus_genes.tsv"),  \
    path("${genome}_summary/${genome}_plasmid_genes.tsv"), emit: genomad_gene_files
    path("${genome}_*")

    script:
    """
    #! /usr/bin/env bash
    set -e
    
    genomad end-to-end \
    --threads ${task.cpus} \
    --cleanup  \
    --enable-score-calibration \
    ${fna} '.' ${params.genomad_db_path} 

    """
}

process RunMEFinder {
    
    label 'mem_medium'
    publishDir "${params.db_path}/pangenomes_annotation/01_mge/${species}/${genome}/mefinder_output", mode: "copy"
    
    input:
    tuple val(genome), val(species), path(fna)

    output:
    tuple val(genome), val(species), path("mefinder.csv"), emit: mefinder_output
    path('*')

    script:
    """
    #! /usr/bin/env bash
    set -e
    export PATH=\$PATH:${params.blastn_path}
    source ${params.resfinder_env_path}/bin/activate
    
    mefinder find --contig ${fna} -t ${task.cpus} mefinder

    """
}

process RunResFinder {
    
    label 'mem_medium'
    publishDir "${params.db_path}/pangenomes_annotation/01_mge/${species}/${genome}/resfinder_output", mode: "copy"
    
    input:
    tuple val(genome), val(species), path(fna)

    output:
    tuple val(genome), val(species), path('ResFinder_results_tab.txt'), emit: resfinder_output
    path('*')

    script:
    """
    #! /usr/bin/env bash
    set -e
    export GIT_PYTHON_GIT_EXECUTABLE=${params.git_executable_path}
    export PATH=\$PATH:${params.git_path}:${params.blastn_path}
    source ${params.resfinder_env_path}/bin/activate

    python -m resfinder \
        -ifa ${fna} \
        -o . \
        -s Other \
        -l ${params.resfinder_min_cov} \
        -t ${params.resfinder_identity_threshold} \
        --acquired \
        -db_res ${params.resfinder_db_path}/resfinder_db \
        -d -db_disinf ${params.resfinder_db_path}/disinfinder_db 

    """
}



process ParsePangenomeAnnotations {

    label 'single_cpu'
    publishDir "${params.db_path}/pangenomes_annotation/03_processed/${species}/", mode: "copy"

    input:
    tuple val(genome), \
    val(species), \
    path(genes_file), \
    path(mefinder_csv), \
    path(resfinder_txt), \
    path(genomad_virus), \
    path(genomad_plasmid), \
    path(contig_len_file), \
    path(eggnog_results), \
    path(max_cluster_info)

    output:
    tuple val(species), \
    path("genomad_virus/${genome}"), \
    path("genomad_plasmid/${genome}"), \
    path("mefinder/${genome}"), \
    path("resfinder/${genome}"), \
    path("eggnog/${genome}")

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    echo ""
    mkdir genomad_virus genomad_plasmid mefinder resfinder eggnog

    python3 ${params.bin_path}/parse_pangenome_annotations.py \
    --species ${species} \
    --genome ${genome} \
    --genes_file ${genes_file} \
    --mefinder_file ${mefinder_csv} \
    --resfinder_file ${resfinder_txt} \
    --genomad_virus_file ${genomad_virus} \
    --genomad_plasmid_file ${genomad_plasmid} \
    --contig_len_file ${contig_len_file} \
    --eggnog_results_file ${eggnog_results} \
    --max_cluster_info_file ${max_cluster_info} 

    """
}

process CombinePangenomeAnnotations {

    label 'single_cpu'
    publishDir "${params.db_path}/pangenomes/${species}/annotation", mode: "copy"

    input:
    tuple val(species), \
    path(genomad_virus_files, stageAs: 'genomad_virus/*'), \
    path(genomad_plasmid_files, stageAs: 'genomad_plasmid/*'), \
    path(mefinder_files, stageAs: 'mefinder/*'), \
    path(resfinder_files, stageAs: 'resfinder/*'),  \
    path(eggnog_files, stageAs: 'eggnog/*')

    output:
    tuple val(species), \
    path('genomad_virus.tsv'), \
    path('genomad_plasmid.tsv'), \
    path('mefinder.tsv'), \
    path('resfinder.tsv'), \
    path('eggnog.tsv')

    shell:
    """
    #! /usr/bin/env bash
    set -e

    COLS_GENOMAD=('gene_id' 'contig_id' 'start' 'end' 'strand' \
    'gene_type' 'contig_length' 'start_genomad' 'end_genomad' \
    'gene' 'annotation_conjscan' 'annotation_amr' \
    'annotation_accessions' 'annotation_description')

    COLS_MEFINDER=('gene_id' 'contig_id' 'start' 'end' 'strand' \
        'gene_type' 'contig_length' 'start_mefinder' 'end_mefinder' \
        'mge_no' 'prediction' 'name' 'type' 'synonyms')

    COLS_RESFINDER=('gene_id' 'contig_id' 'start' 'end' 'strand' \
        'gene_type' 'contig_length' 'start_resfinder' 'end_resfinder' \
        'resistance_gene' 'phenotype' 'accession_no')

    COLS_EGGNOG=('#query' 'seed_ortholog' 'evalue' 'score' 'eggNOG_OGs' \
        'max_annot_lvl' 'COG_category' 'Description' 'Preferred_name' 'GOs' \
        'EC' 'KEGG_ko' 'KEGG_Pathway' 'KEGG_Module' 'KEGG_Reaction' \
        'KEGG_rclass' 'BRITE' 'KEGG_TC' 'CAZy' 'BiGG_Reaction' 'PFAMs' \
        'contig_id' 'start' 'end' 'strand' 'gene_type' 'contig_length')

    # Join array elements with tabs and append a newline
    echo -e "\$(IFS=\$'\\t'; echo \"\${COLS_GENOMAD[*]}\")" > genomad_virus.tsv
    cat !{genomad_virus_files} >> genomad_virus.tsv

    echo -e "\$(IFS=\$'\\t'; echo \"\${COLS_GENOMAD[*]}\")" > genomad_plasmid.tsv
    cat !{genomad_plasmid_files} >> genomad_plasmid.tsv

    echo -e "\$(IFS=\$'\\t'; echo \"\${COLS_MEFINDER[*]}\")" > mefinder.tsv
    cat !{mefinder_files} >> mefinder.tsv

    echo -e "\$(IFS=\$'\\t'; echo \"\${COLS_RESFINDER[*]}\")" > resfinder.tsv
    cat !{resfinder_files} >> resfinder.tsv

    echo -e "\$(IFS=\$'\\t'; echo \"\${COLS_EGGNOG[*]}\")" > eggnog.tsv
    cat !{eggnog_files} >> eggnog.tsv

    """
}


process EnhancePangenome {

    label 'single_cpu'
    publishDir (
        path: "${params.db_path}/pangenomes/${species}",
        mode: "copy",
        pattern: "*.tsv",
        saveAs: { fn ->
            if (fn.endsWith("annot.tsv")) { "annotation/${fn}" }
            else { "${fn}" }
        }
    )

    input:
    tuple val(species), \
    path(genomad_virus), \
    path(genomad_plasmid), \
    path(mefinder_file), \
    path(resfinder_file), \
    path(eggnog_file), \
    path(genes_info_file), \
    path(cluster_info_file_list, stageAs: 'old/*')

    output:
    tuple val(species), \
    path("genes_annotated.tsv"), \
    path('clusters_*_annot.tsv'), \
    path('clusters_*_info.tsv')

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    cluster_pct_list=\$(echo ${params.centroid_cluster_percents} | tr -d '[],')

    python3 ${params.bin_path}/enhance_pangenome.py \
    --cluster_thresholds \${cluster_pct_list}  \
    --genomad_virus_file ${genomad_virus} \
    --genomad_plasmid_file ${genomad_plasmid} \
    --mefinder_file ${mefinder_file} \
    --resfinder_file ${resfinder_file} \
    --eggnog_file ${eggnog_file} \
    --genes_info_file ${genes_info_file} \
    --cluster_info_file_list ${cluster_info_file_list} 

    """
}

process ComputeRunSNPsChunks {
    label 'single_cpu'
    publishDir "${params.db_path}/chunks/sites/run/chunksize.${params.run_chunk_size}/${species}/", mode: "copy"

    input:
    tuple val(genome), val(species), path(fna)

    output:
    path("${genome}.json")

"""
#!/usr/bin/env python3

import json
from midas.models.species import design_run_snps_chunks, design_merge_snps_chunks
from midas.common.utils import OutputStream

run_snp_chunks_to_cache = design_run_snps_chunks('${species}', '${fna}', ${params.run_chunk_size})

with OutputStream('${genome}.json') as stream:
        json.dump(run_snp_chunks_to_cache, stream)

"""
}

process ComputeMergeSNPsChunks {
    label 'single_cpu'
    publishDir "${params.db_path}/chunks/sites/merge/chunksize.${params.merge_chunk_size}/${species}/", mode: "copy"

    input:
    tuple val(genome), val(species), path(fna)

    output:
    path("${genome}.json")

"""
#!/usr/bin/env python3

import json
from midas.models.species import design_merge_snps_chunks
from midas.common.utils import OutputStream

merge_snp_chunks_to_cache = design_merge_snps_chunks('${species}', '${fna}', ${params.merge_chunk_size})

with OutputStream('${genome}.json') as stream:
        json.dump(merge_snp_chunks_to_cache, stream)

"""
}
