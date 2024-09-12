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

workflow {

    genomes = Channel
        .fromPath(params.genomes_tsv_path)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.genome, row.species, row.representative, row.genome_is_representative) } 

    // rep_genomes = Channel
    //     .fromPath(params.genomes_tsv_path)
    //     .splitCsv(header: true, sep: '\t')
    //     .filter { r -> (r.genome_is_representative == "1")}
    //     .map { row -> tuple(row.genome, row.species, row.representative, row.genome_is_representative) } 

    // rep_genomes.view()
    
    annotateGenomes(genomes)
    generateGeneFeatures(annotateGenomes.out.genome, annotateGenomes.out.species, annotateGenomes.out.gff)
    hmmMarkerSearch(annotateGenomes.out.genome, annotateGenomes.out.species, annotateGenomes.out.faa, annotateGenomes.out.ffn)
    inferMarkers(hmmMarkerSearch.out)
    
    buildMarkerDB(inferMarkers.out.markers_fa.collect(), inferMarkers.out.markers_map.collect())
    // x.view()
    
    // buildMarkerDB(repgenomes)
        

}


process annotateGenomes {
    label 'mem_low'
    errorStrategy 'finish'
    conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/mtest'
    publishDir "${params.db_dir_path}/gene_annotations/${species}/${genome}", mode: "copy"

    input:
    tuple val(genome), val(species), val(representative), val(genome_is_representative)

    output:
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
    val(genome)
    val(species)
    path(hmmsearch)
    path(ffn)

    output:
    val(genome), emit: genome
    val(species), emit: species
    path("${genome}.markers.fa"), emit: markers_fa
    path("${genome}.markers.map"), emit: markers_map


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


    // #! /usr/bin/env bash
    // set -e
    // set -x

    // #get rep genome for each species

    
    // #example test_db/markers/phyeco/temp/117086/GCA_900552055.1/GCA_900552055.1.markers.fa
    // for f in *.fa; 
    // do cat "$f" >> phyeco.fa; done



    // for f in *.map; do cat "$f" >> phyeco.map; done

    // hs-blastn index phyeco.fa 

// process buildPangenomeCommand {
//     label 'mem_medium'
//     errorStrategy 'finish'
//     conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/MIDASv3'

//     input:
//     tuple val(genome), file from build_markerdb_out, val(db_name), file(db_dir)

//     output:
//     tuple val(genome), file("build_pangenome_output/${genome}") into build_pangenome_out

//     script:
//     """
//     #! /usr/bin/env bash
//     set -e
//     set -x

//     conda activate MIDASv3

//       scratch_dir="."

//     midas build_pangenome --species all --midasdb_name ${db_name} --midasdb_dir ${db_dir} --debug --force --scratch_dir \${scratch_dir} -t ${task.cpus} --recluster
//     """
// }

// process pangenomePerSpecies {
//     label 'mem_medium'
//     errorStrategy 'finish'
//     conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/MIDASv3'

//     input:
//     tuple val(genome), file from build_pangenome_out, val(db_name), file(db_dir)

//     output:
//     tuple val(genome), file("pangenome_per_species_output/${genome}") into pangenome_per_species_out

//     script:
//     """
//     #! /usr/bin/env bash
//     set -e
//     set -x

//     conda activate MIDASv3

//       species_list="\$db_dir/species_list.txt"

//     species=\$(awk "NR==\$SGE_TASK_ID" \$species_list | awk -F',' '{print \$1}')
//     echo \$species

//     bash /wynton/protected/home/sirota/clairedubin/bin/MIDAS/bin/pipeline.sh \$species \$db_dir/pangenomes \$${task.cpus} 8000 /wynton/protected/home/sirota/clairedubin/bin/MIDAS/bin
//     """
// }

// process reclusterCentroids {
//     label 'mem_medium'
//     errorStrategy 'finish'
//     conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/MIDASv3'

//     input:
//     tuple val(genome), file from pangenome_per_species_out, val(db_name), file(db_dir)

//     output:
//     tuple val(genome), file("recluster_centroids_output/${genome}") into recluster_centroids_out

//     script:
//     """
//     #! /usr/bin/env bash
//     set -e
//     set -x

//     conda activate MIDASv3

//       scratch_dir="./scratch"

//     midas recluster_centroids --species all --midasdb_name ${db_name} --midasdb_dir ${db_dir} --debug --force --scratch_dir \${scratch_dir} -t ${task.cpus}
//     """
// }

// process augmentPangenome {
//     label 'mem_medium'
//     errorStrategy 'finish'
//     conda '/wynton/protected/home/sirota/clairedubin/anaconda3/envs/MIDASv3'

//     input:
//     tuple val(genome), file from recluster_centroids_out, val(db_name), file(db_dir)

//     output:
//     tuple val(genome), file("augment_pangenome_output/${genome}") into augment_pangenome_out

//     script:
//     """
//     #! /usr/bin/env bash
//     set -e
//     set -x

//     conda activate MIDASv3

//     scratch_dir="./scratch"

//     midas augment_pangenome --species all --midasdb_name ${db_name} --midasdb_dir ${db_dir} --debug --force --scratch_dir \${scratch_dir} -t ${task.cpus}
//     """
// }

// process annotatePerGenome {
//     label 'mem_medium'
//     errorStrategy 'finish'

//     input:
//     tuple val(genome), file from augment_pangenome_out, val(db_name), file(db_dir)

//     output:
//     tuple val(genome), file("annotate_per_genome_output/${genome}") into annotate_per_genome_out

//     script:
//     """
//     #! /usr/bin/env bash
//     set -e

//     midasdb_dir="${db_dir}"
//     genomes_file="\${midasdb_dir}/genomes.tsv"
//     genome_id=\$(awk "NR==\$SGE_TASK_ID" \$genomes_file | awk -F'\\t' '{print \$1}')
//     species_id=\$(awk "NR==\$SGE_TASK_ID" \$genomes_file | awk -F'\\t' '{print \$2}')

//     GENOMAD_DATA_DIR="/wynton/protected/home/sirota/clairedubin/databases/genomad_db"
//     RESFIND_DATA_DIR="/wynton/protected/home/sirota/clairedubin/databases/resfinder_dbs"

//     echo "Running GENOMAD: \$species_id:\$genome_id"

//     conda activate genomad

//     ffn="\${m
