process ClusterCentroids {
    label 'mem_high'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/mtest'
    publishDir "${params.db_dir_path}/pangenomes/${species}/temp/vsearch/", mode: "copy"

    input:
    tuple val(cluster_pct), val(species), path(genes_ffn)

    output:
    tuple val(species), path("centroids.${cluster_pct}.ffn"), emit: centroid_ffn
    tuple val(species), path("uclust.${cluster_pct}.txt"), emit: uclust
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

process ReClusterCentroids {
    label 'mem_high'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/mtest'
    publishDir "${params.db_dir_path}/pangenomes/${species}/temp/cdhit", mode: "copy"

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
