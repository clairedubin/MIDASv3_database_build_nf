#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.max_cluster_val = params.centroid_cluster_percents.max()
params.rarefaction_rep_counts = 10

params.original_db_path = "${params.db_output_dir}/${params.db_name}/"
params.rarefaction_output_path = "${params.db_output_dir}/${params.db_name}/rarefaction/"
// New optional parameter, default to null
params.species_sample_size_csv = null 

vsearch_container = "quay.io/biocontainers/vsearch:2.30.0--hd6d6fdc_0"
pandas_container = "quay.io/biocontainers/pandas:2.2.1"
biopython_container = "quay.io/biocontainers/biopython:1.70--np112py27_1"

workflow {

    genomes = Channel
        .fromPath(params.genomes_tsv_path)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.genome, row.species, row.representative, row.genome_is_representative, row.fasta_path) } 

    genomes
        .map { r -> r[1] }  
        .unique()           
        .combine(Channel.of(params.genomes_tsv_path))
        .set { species_list }

    // --- NEW LOGIC START ---
    // Prepare inputs for MakeRarefactionInputs based on whether optional CSV exists
    if (params.species_sample_size_csv) {
        Channel
            .fromPath(params.species_sample_size_csv)
            .splitCsv(header: true)
            .map { row -> tuple(row.species, row.starting_sample_size) }
            .set { size_ch }

        // Join species list with sizes. 
        // remainder: true ensures species NOT in the CSV still get passed through with a null size
        species_list
            .join(size_ch, remainder: true)
            .set { make_rarefaction_input_ch }
    } else {
        // If no CSV provided, pass null as the third element for all species
        species_list
            .map { species, tsv -> tuple(species, tsv, null) }
            .set { make_rarefaction_input_ch }
    }
    // --- NEW LOGIC END ---

    MakeRarefactionInputs(make_rarefaction_input_ch)

    def rarefaction_inputs = MakeRarefactionInputs.out.transpose()
        .map { species_code, csv_path ->
        def filename = csv_path.toString().tokenize('/')[-1]  // Extract filename from path
        def match = filename =~ /${species_code}_(\d+)genomes_rep(\d+)\.csv/  // Regex to extract numbers
        if (!match) {
            error "Filename format unexpected: ${filename}"
        }
        def (species_num, rep_num) = match[0][1..2]

        return tuple(species_code, species_num.toInteger(), rep_num.toInteger(), csv_path)
    }

    CombineCleanedGenes(rarefaction_inputs)

    Channel.of(params.max_cluster_val)
        .combine(CombineCleanedGenes.out.genes_ffn)
        .set{max_cluster_ch}
    
    ClusterCentroids(max_cluster_ch)
        .set {max_cluster_output}

    CleanCentroids(
        max_cluster_output.centroid_ffn,
        max_cluster_output.cluster_pct
    )

    // Clustering results of highest cluster threshold (e.g. C99) at lower thresholds
    remaining_clusters_list = params.centroid_cluster_percents.findAll {it != params.max_cluster_val}
        
    Channel
        .fromList(remaining_clusters_list)
        .combine(max_cluster_output.centroid_ffn)
        .set{lower_cluster_ch}

    ClusterCentroidsLowerThresholds(lower_cluster_ch)
        .set {lower_cluster_output}

    lower_cluster_output.uclust
        .concat(max_cluster_output.uclust)
        .groupTuple(by: [0,1,2])
        .set {uclust_to_parse}

    ParseCentroidInfo(uclust_to_parse)

    ParseCentroidInfo.out
        .join(CleanCentroids.out.cleaned_centroids, by: [0,1,2])
        .join(CombineCleanedGenes.out.genes_ffn, by: [0,1,2])
        .join(CombineCleanedGenes.out.genes_len, by: [0,1,2])
        .set{refine_clusters_input}

    RefineClusters(refine_clusters_input)

    // Redoing clustering at lower thresholds        
    Channel
        .fromList(remaining_clusters_list)
        .combine(RefineClusters.out.centroid_ffn)
        .set{lower_cluster_refine_ch}

    ReClusterCentroids(lower_cluster_refine_ch)

    RefineClusters.out.gene_files
        .combine(ReClusterCentroids.out.uclust.groupTuple(by: [0,1,2]), by: [0,1,2])
        .set{recluster_to_parse}

    ParseReclusteredCentroidInfo(recluster_to_parse)

    CombinePangenomeSizes(ParseReclusteredCentroidInfo.out.pangenome_size.collect())

}

process MakeRarefactionInputs {

    label 'single_cpu'
    publishDir "${params.rarefaction_output_path}/${species}/subsamples/", mode: "copy"
    errorStrategy 'terminate'

    input:
    // Updated input tuple to accept optional max_size
    tuple val(species), path(genomes_tsv), val(max_size)

    output:
    tuple val(species), path("${species}_*genomes_rep*.csv")

    script:
    // Groovy logic to conditionally create the argument string
    def max_size_arg = max_size ? "--max_size ${max_size}" : ""
    
    """
    #! /usr/bin/env bash
    set -e
    set -x

   make_rarefaction_input.py \
   --genomes_tsv ${genomes_tsv} \
   --species ${species} \
   --num_reps ${params.rarefaction_rep_counts} \
   --db_path ${params.original_db_path} \
   ${max_size_arg}
    """

}

process CombineCleanedGenes {
    label 'single_cpu'
    errorStrategy 'terminate'
    conda false

    input:
    tuple val(species), val(num_genomes), val(rep_num), path(rarefaction_csv)

    output:
    tuple val(species), val(num_genomes), val(rep_num), path("genes.ffn"), emit: genes_ffn
    tuple val(species), val(num_genomes), val(rep_num), path("genes.len"), emit: genes_len

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    gene_ffns=\$(awk -F',' 'NR>1 {print \$3}' ${rarefaction_csv})
    gene_lens=\$(awk -F',' 'NR>1 {print \$4}' ${rarefaction_csv})

    cat \$gene_ffns > genes.ffn
    cat \$gene_lens > genes.len
    
    """
}

process ClusterCentroids {
    label 'mem_high'
    container "${vsearch_container}"
    conda false

    input:
    tuple val(cluster_pct), val(species), val(num_genomes), val(rep_num), path(genes_ffn)

    output:
    tuple val(species), val(num_genomes), val(rep_num), path("centroids.${cluster_pct}.ffn"), emit: centroid_ffn
    tuple val(species), val(num_genomes), val(rep_num), path("uclust.${cluster_pct}.txt"), emit: uclust
    val cluster_pct, emit: cluster_pct


    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    if [ ! -s ${genes_ffn} ]; then
        echo "ERROR: ${genes_ffn} is empty"
        exit 1 
    fi

    cluster_prop="\$(awk "BEGIN {printf \\"%.2f\\", ${cluster_pct} / 100}")"
    vsearch --cluster_fast ${genes_ffn} --id \${cluster_prop} --threads ${task.cpus} --centroids centroids.${cluster_pct}.ffn --uc uclust.${cluster_pct}.txt

    """
}

process ClusterCentroidsLowerThresholds {
    label 'mem_high'
    container "${vsearch_container}"
    conda false

    input:
    tuple val(cluster_pct), val(species), val(num_genomes), val(rep_num), path(genes_ffn)

    output:
    tuple val(species), val(num_genomes), val(rep_num), path("centroids.${cluster_pct}.ffn"), emit: centroid_ffn
    tuple val(species), val(num_genomes), val(rep_num), path("uclust.${cluster_pct}.txt"), emit: uclust
    val cluster_pct, emit: cluster_pct


    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    if [ ! -s ${genes_ffn} ]; then
        echo "ERROR: ${genes_ffn} is empty"
        exit 1 
    fi

    cluster_prop="\$(awk "BEGIN {printf \\"%.2f\\", ${cluster_pct} / 100}")"
    vsearch --cluster_fast ${genes_ffn} --id \${cluster_prop} --threads ${task.cpus} --centroids centroids.${cluster_pct}.ffn --uc uclust.${cluster_pct}.txt

    """
}

process CleanCentroids {
    label 'single_cpu'
    errorStrategy 'terminate'

    input:
    tuple val(species), val(num_genomes), val(rep_num), path(centroid_ffn)
    val(cluster_pct)

    output:
    tuple val(species), val(num_genomes), val(rep_num), \
    path("centroids.${cluster_pct}.clean.ffn"), \
    path("centroids.${cluster_pct}.ambiguous.ffn"), emit: cleaned_centroids

"""
#!/usr/bin/env python3

import Bio.SeqIO
import os
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
    errorStrategy 'retry'

    input:
    tuple val(species), val(num_genomes), val(rep_num), path(uclust_files)

    output:
    tuple val(species), val(num_genomes), val(rep_num), path("gene_info.txt")

    script:
    """
    parse_centroids.py ${uclust_files}
    """
}

process RefineClusters {

    label 'mem_medium'
    errorStrategy 'retry'

    
    input:
    tuple val(species), val(num_genomes), val(rep_num), path(gene_info), path(clean_centroids), path(amb_centroids), path(genes_ffn), path(genes_len)

    output:
    tuple val(species), val(num_genomes), val(rep_num), path("centroids.ffn"), emit: centroid_ffn 
    tuple val(species), val(num_genomes), val(rep_num), path("temp/cdhit/genes.len"), path("temp/cdhit/gene_info.txt"), emit: gene_files
    path("temp/cdhit/*")
    path("temp/cdhit/centroids.${params.max_cluster_val}.ffn")

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x
    
    pipeline.sh \
        ${species} \
        ${gene_info} \
        ${clean_centroids} \
        ${amb_centroids} \
        ${genes_ffn} \
        ${genes_len} \
        ${task.cpus} \
        ${task.memory.toMega()} 

    cp "temp/cdhit/centroids.${params.max_cluster_val}.ffn" centroids.ffn
    """
}

process ReClusterCentroids {
    label 'mem_high'
    errorStrategy 'retry'
    container "${vsearch_container}"
    conda false
    
    input:
    tuple val(cluster_pct), val(species), val(num_genomes), val(rep_num), path(genes_ffn)


    output:
    tuple val(species), val(num_genomes), val(rep_num), path("centroids.${cluster_pct}.ffn"), emit: centroid_ffn
    tuple val(species), val(num_genomes), val(rep_num), path("uclust.${cluster_pct}.txt"), emit: uclust


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
    publishDir "${params.rarefaction_output_path}/${species}/${num_genomes}/${rep_num}/", mode: "copy"
    errorStrategy 'retry'

    input:
    tuple val(species), val(num_genomes), val(rep_num), path(genes_len), path(genes_info), path(uclust_files)

    output:
    tuple val(species), val(num_genomes), val(rep_num), path('genes_info.tsv'), emit: reclustered_centroid_info
    path("pangenome_size_${species}_${num_genomes}_${rep_num}.csv"), emit: pangenome_size


    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    mkdir temp
    cat ${params.original_db_path}/markers/${params.marker_set}/temp/${species}/*/*.markers.map >> temp/mapfile

    parse_reclustered_centroids.py \
    --gene_info_file ${genes_info} \
    --max_percent ${params.max_cluster_val} \
    --gene_length_file ${genes_len} \
    --uclust_files ${uclust_files} 

    
    calc_pangenome_size.py \
    --gene_info_file genes_info.tsv \
    --species ${species} \
    --num_genomes ${num_genomes} \
    --rep_number ${rep_num}
    """
}

process CombinePangenomeSizes {

    label 'single_cpu'
    publishDir "${params.rarefaction_output_path}/", mode: "copy"
    errorStrategy 'terminate'
    conda false

    input:
    path('pangenome_size*.tsv')

    output:
    path('rarified_pangenome_sizes.csv')
    
    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    { head -n 1 -q pangenome_size*.tsv | head -n 1; tail -n +2 -q pangenome_size*.tsv; } > rarified_pangenome_sizes.csv
    """
}
