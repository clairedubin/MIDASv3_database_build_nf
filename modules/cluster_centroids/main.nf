process ClusterCentroids {
    label 'mem_medium'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/mtest'
    publishDir "${params.db_dir_path}/pangenomes/${species}/temp/vsearch/"

    input:
    tuple val(cluster_pct), val(species), path(genes_ffn), path(genes_len)


    output:
    tuple val(species), path("centroids.${cluster_pct}.ffn"), path(genes_len), emit: recluster_input
    val species, emit: species
    val cluster_pct, emit: cluster_pct
    path("centroids.${cluster_pct}.ffn"), emit: centroids_ffn
    path(genes_len), emit: genes_len
    path("uclust.${cluster_pct}.txt"), emit: uclust



    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x

    cluster_prop="\$(awk "BEGIN {printf \\"%.2f\\", ${cluster_pct} / 100}")"
    vsearch --cluster_fast ${genes_ffn} --id \${cluster_prop} --threads ${task.cpus} --centroids centroids.${cluster_pct}.ffn --uc uclust.${cluster_pct}.txt


    """
}