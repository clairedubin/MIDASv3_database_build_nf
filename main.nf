#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- Parameters & Path Setup ---
params.max_cluster_val = params.centroid_cluster_percents.max()

// Ensure directories end with trailing "/" characters and standardize paths
params.findAll { key, _ -> key.endsWith("_dir")}
    .each { key, path ->
        if (!path.endsWith("/")) {
            params[key.replace("_dir", "_path")] = "${path}/"
        } else {
            params[key.replace("_dir", "_path")] = "${path}"
        }
    }

params.db_path = "${params.db_output_dir}/${params.db_name}/"
params.marker_model_hmm_path = "${baseDir}/bin/markers_models/${params.marker_set}/marker_genes.hmm"

// Ensure the main output directory exists
def output_db = file(params.db_path)
if (!output_db.exists()) {
    output_db.mkdirs()
}

// ---------------------------------------------------------
// INCLUDE MODULES
// ---------------------------------------------------------
include { check_input } from './modules/check_input'
include { ClusterCentroids as ClusterCentroids } from './modules/cluster_centroids'
include { ClusterCentroids as ClusterCentroidsLowerThresholds } from './modules/cluster_centroids'

workflow {

    //// 1. PREPARE DATABASES (Downloads only if missing)
    DownloadGeNomadDB()
    DownloadEggNogDB()
    DownloadResFinderDB()

    //// 2. CHECK INPUT FILE FORMATTING
    check_input(params.genomes_tsv_path)

    genomes = Channel
        .fromPath(params.genomes_tsv_path)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.genome, row.species, row.representative, row.genome_is_representative, row.fasta_path) }

    // Filter to only representative genomes
    genomes.filter { r -> (r[3] == "1") }
        .map{ r -> tuple(r[0], r[1]) }
        .set{ rep_genomes }

    //// 3. IDENTIFY AND ANNOTATE GENES
    AnnotateGenomes(genomes)
    CombineCleanedGenes(AnnotateGenomes.out.cleaned_genes.groupTuple(by: 1))
    CalculateContigLength(
        AnnotateGenomes.out.fna_tuple
            .groupTuple(by: 1)
            .map{ r -> tuple(r[1], r[2])}
    )

    //// 4. COMPUTE CHUNKS
    rep_genomes
        .join(AnnotateGenomes.out.fna_tuple, by: [0,1])
        .set{rep_genome_fna_tuples}
    ComputeChunks(rep_genome_fna_tuples)

    //// 5. BUILD MARKER DB
    HMMMarkerSearch(AnnotateGenomes.out.faa_ffn_tuple)
    rep_genomes
        .join(HMMMarkerSearch.out.markers, by: [0,1])
        .set{ rep_inferred_markers }
    BuildMarkerDB(
        rep_inferred_markers.collect{ r -> r[2]},
        rep_inferred_markers.collect{ r -> r[3]}
    )

    //// 6. CLUSTERING
    Channel.of(params.max_cluster_val)
        .combine(CombineCleanedGenes.out.genes_ffn)
        .set{max_cluster_ch}
    
    ClusterCentroids(max_cluster_ch)
        .set {max_cluster_output}

    CleanCentroids(
        max_cluster_output.centroid_ffn,
        max_cluster_output.cluster_pct
    )

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

    Channel
        .fromList(remaining_clusters_list)
        .combine(RefineClusters.out.centroid_ffn)
        .set{lower_cluster_refine_ch}
    
    ReClusterCentroids(lower_cluster_refine_ch)

    HMMMarkerSearch.out.markers
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

    //// 7. FUNCTIONAL ANNOTATION (Linked to Corrected DB Paths)
    RunEggNog(RefineClusters.out.centroid_ffn, DownloadEggNogDB.out.db)
    RunGeNomad(AnnotateGenomes.out.fna_tuple, DownloadGeNomadDB.out.db)
    
    RunResFinder(
        AnnotateGenomes.out.fna_tuple, 
        DownloadResFinderDB.out.res_db, 
        DownloadResFinderDB.out.dis_db
    )

    //// 8. COMBINE RESULTS
    AnnotateGenomes.out.genes
        .combine(RunResFinder.out.resfinder_output, by: [0,1])
        .combine(RunGeNomad.out.genomad_gene_files, by: [0,1])
        .map{it.swap(1,0)}
        .combine(CalculateContigLength.out, by: 0)
        .combine(RunEggNog.out.eggnog_output, by: 0)
        .combine(AugmentPangenomes.out.max_cluster_info, by: 0)
        .map{it.swap(1,0)}
        .set{annotations_to_parse}
    
    ParsePangenomeAnnotations(annotations_to_parse)
    CombinePangenomeAnnotations(ParsePangenomeAnnotations.out.groupTuple(by: 0))

    CombinePangenomeAnnotations.out
        .combine(ParseReclusteredCentroidInfo.out.reclustered_centroid_info, by: 0)
        .combine(AugmentPangenomes.out.cluster_info_file_list, by: 0)
        .set{enhance_pangenome_input}
    
    EnhancePangenome(enhance_pangenome_input)
}

// ─────────────────────────────────────────────────────────────────────────────
// DATABASE SETUP PROCESSES (mkdir fix applied)
// ─────────────────────────────────────────────────────────────────────────────


process DownloadGeNomadDB {
    executor 'local'
    conda false
    storeDir "${params.genomad_db_dir}"

    output:
    path "genomad_db", emit: db

    script:
    """
    # Download the database archive using the URL from params
    wget -O genomad_db_*.tar.gz "${params.genomad_db_url}"
    
    # Extract the archive
    tar -xzf genomad_db_*.tar.gz
    
    # Clean up the tarball
    rm genomad_db_*.tar.gz
    """
}

process DownloadEggNogDB {
    executor 'local'
    conda false
    storeDir "${params.eggnog_db_dir}"

    output:
    path "data", emit: db

    script:
    def base_url = "http://eggnog6.embl.de/download/emapperdb-5.0.2"
    """
    mkdir -p data && cd data
    wget --no-cache --show-progress ${base_url}/eggnog.db.gz && gunzip eggnog.db.gz
    wget --no-cache --show-progress ${base_url}/eggnog_proteins.dmnd.gz && gunzip eggnog_proteins.dmnd.gz
    wget --no-cache --show-progress ${base_url}/eggnog.taxa.tar.gz && tar -xzvf eggnog.taxa.tar.gz && rm eggnog.taxa.tar.gz
    wget --no-cache --show-progress ${base_url}/pfam.tar.gz && mkdir -p pfam && tar -xzvf pfam.tar.gz -C pfam --strip-components=1 && rm pfam.tar.gz
    """
}

process DownloadResFinderDB {
    executor 'local'
    conda "${projectDir}/env/resfinder.yml"
    storeDir "${params.resfinder_db_dir}"

    output:
    path "resfinder_db",   emit: res_db
    path "pointfinder_db", emit: point_db
    path "disinfinder_db", emit: dis_db

    script:
    """
    git clone -b 2.4.0 https://git@bitbucket.org/genomicepidemiology/resfinder_db.git resfinder_db    
    git clone -b 4.1.1 https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git pointfinder_db
    git clone -b 2.0.1 https://git@bitbucket.org/genomicepidemiology/disinfinder_db.git disinfinder_db
    
    cd resfinder_db && python3 INSTALL.py
    cd ../pointfinder_db && python3 INSTALL.py
    cd ../disinfinder_db && python3 INSTALL.py
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// CORE PIPELINE PROCESSES
// ─────────────────────────────────────────────────────────────────────────────

process AnnotateGenomes {
    label 'mem_medium'
    conda "${projectDir}/env/prokka.yml"
    publishDir "${params.db_path}/gene_annotations/${species}/${genome}", mode: "copy"

    input:
    tuple val(genome), val(species), val(representative), val(genome_is_representative), val(genome_path)

    output:
    tuple val(genome), val(species), path("${genome}.fna"),                         emit: fna_tuple
    tuple val(genome), val(species), path("${genome}.faa"), path("${genome}.ffn"), emit: faa_ffn_tuple
    tuple val(genome), val(species), path("${genome}.genes"),                      emit: genes
    tuple val(genome), val(species), path("${genome}.genes.ffn"), path("${genome}.genes.len"), emit: cleaned_genes
    path("${genome}.gff.db")
    path("${genome}.tsv")
    path("${genome}.gff")

    script:
    """
    prokka --kingdom Bacteria --metagenome --cpus ${task.cpus} --prefix ${genome} --locustag ${genome} --outdir . --compliant --force ${genome_path}
    generate_gene_features.py --gff ${genome}.gff --genome ${genome}
    clean_genes.py --ffn ${genome}.ffn --genome ${genome}
    """
}

process HMMMarkerSearch {
    label 'mem_medium'
    conda "${projectDir}/env/hmmer.yml"
    publishDir "${params.db_path}/markers/${params.marker_set}/temp/${species}/${genome}", mode: "copy"

    input:
    tuple val(genome), val(species), path(faa), path(ffn)

    output:
    tuple val(genome), val(species), path("${genome}.markers.fa"), path("${genome}.markers.map"), emit: markers
    path("${genome}.hmmsearch")

    script:
    """
    hmmsearch --noali --cpu ${task.cpus} --domtblout "${genome}.hmmsearch" ${params.marker_model_hmm_path} ${faa}
    infer_markers.py --genome ${genome} --species ${species} --hmmsearch_file ${genome}.hmmsearch --annotation_ffn ${ffn} --hmmsearch_max_evalue ${params.hmmsearch_max_evalue} --hmmsearch_min_cov ${params.hmmsearch_min_cov}
    """
}

process BuildMarkerDB {
    label 'single_cpu'
    conda "${projectDir}/env/hsblastn.yml"
    publishDir "${params.db_path}/markers/${params.marker_set}", mode: "copy"

    input:
    path '*.markers.fa'
    path '*.markers.map'

    output:
    path("${params.marker_set}.*")

    script:
    """
    cat *.markers.fa > ${params.marker_set}.fa
    cat *.markers.map > ${params.marker_set}.map
    hs-blastn index ${params.marker_set}.fa
    """
}

process CombineCleanedGenes {
    label 'single_cpu'
    publishDir "${params.db_path}/pangenomes/${species}", mode: "copy"

    input:
    tuple val(genome), val(species), path(gene_ffns), path(gene_lens)

    output:
    tuple val(species), path("genes.ffn"), emit: genes_ffn
    tuple val(species), path("genes.len"), emit: genes_len
    tuple val(species), path("temp/vsearch/genes.len")

    script:
    """
    cat ${gene_ffns} > genes.ffn
    cat ${gene_lens} > genes.len
    mkdir -p temp/vsearch
    cp genes.len temp/vsearch/genes.len
    """
}

process CleanCentroids {
    label 'single_cpu'
    conda "${projectDir}/env/python-bio.yml"
    publishDir "${params.db_path}/pangenomes/${species}/temp/vsearch/", mode: "copy"

    input:
    tuple val(species), path(centroid_ffn)
    val(cluster_pct)

    output:
    tuple val(species), \
        path("centroids.${cluster_pct}.clean.ffn"), \
        path("centroids.${cluster_pct}.ambiguous.ffn"), emit: cleaned_centroids

    script:
    """
    #!/usr/bin/env python3
    import Bio.SeqIO
    import os
    from common import has_ambiguous_bases

    ffn_in = "${centroid_ffn}"
    assert os.path.getsize(ffn_in) > 0

    output_ambiguous = ffn_in.replace('.ffn', '.ambiguous.ffn')
    output_clean     = ffn_in.replace('.ffn', '.clean.ffn')

    with open(output_ambiguous, 'w') as o_ambiguous, \\
         open(output_clean, 'w') as o_clean:
        for rec in Bio.SeqIO.parse(ffn_in, 'fasta'):
            c_id  = rec.id
            c_seq = str(rec.seq).upper()
            if has_ambiguous_bases(c_seq):
                o_ambiguous.write(f">{c_id}\\n{c_seq}\\n")
            else:
                o_clean.write(f">{c_id}\\n{c_seq}\\n")
    """
}

process ParseCentroidInfo {
    label 'single_cpu'
    conda "${projectDir}/env/python-bio.yml"
    publishDir "${params.db_path}/pangenomes/${species}/temp/vsearch/", mode: "copy"

    input:
    tuple val(species), path(uclust_files)

    output:
    tuple val(species), path("gene_info.txt")

    script:
    """
    parse_centroids.py ${uclust_files}
    """
}

process RefineClusters {
    label 'mem_medium'
    conda "${projectDir}/env/cdhit.yml"
    publishDir ( path: "${params.db_path}/pangenomes/${species}", mode: "copy" )

    input:
    tuple val(species), path(gene_info), path(clean_centroids), path(amb_centroids), path(genes_ffn), path(genes_len)

    output:
    tuple val(species), path("centroids.ffn"),                                        emit: centroid_ffn
    tuple val(species), path("temp/cdhit/genes.len"), path("temp/cdhit/gene_info.txt"), emit: gene_files
    path("temp/cdhit/*")

    script:
    """
    pipeline.sh ${species} ${gene_info} ${clean_centroids} ${amb_centroids} ${genes_ffn} ${genes_len} ${task.cpus} ${task.memory.toMega()}
    cp "temp/cdhit/centroids.${params.max_cluster_val}.ffn" centroids.ffn
    """
}

process ReClusterCentroids {
    label 'mem_high'
    conda "${projectDir}/env/vsearch.yml"
    publishDir "${params.db_path}/pangenomes/${species}/temp/", mode: "copy"

    input:
    tuple val(cluster_pct), val(species), path(genes_ffn)

    output:
    tuple val(species), path("centroids.${cluster_pct}.ffn"), emit: centroid_ffn
    tuple val(species), path("uclust.${cluster_pct}.txt"),    emit: uclust

    script:
    """
    cluster_prop="\$(awk "BEGIN {printf \\"%.2f\\", ${cluster_pct} / 100}")"
    vsearch --cluster_fast ${genes_ffn} --id \${cluster_prop} --threads ${task.cpus} --centroids centroids.${cluster_pct}.ffn --uc uclust.${cluster_pct}.txt
    """
}

process ParseReclusteredCentroidInfo {
    label 'single_cpu'
    conda "${projectDir}/env/python-bio.yml"
    publishDir "${params.db_path}/pangenomes/${species}/", mode: "copy"

    input:
    tuple val(species), path(genes_len), path(genes_info), path(uclust_files), path(marker_map_files)

    output:
    tuple val(species), path('genes_info.tsv'), emit: reclustered_centroid_info
    path('temp/*')

    script:
    """
    mkdir temp
    cat *.markers.map >> temp/mapfile
    parse_reclustered_centroids.py --gene_info_file ${genes_info} --max_percent ${params.max_cluster_val} --gene_length_file ${genes_len} --uclust_files ${uclust_files}
    cp genes_info.tsv temp/genes_info.tsv
    cp temp/mapfile temp/markers.map
    """
}

process AugmentPangenomes {
    label 'single_cpu'
    conda "${projectDir}/env/python-bio.yml"
    publishDir "${params.db_path}/pangenomes/${species}/augment/", mode: "copy"

    input:
    tuple val(species), val(rep_genome), path(gene_info_tsv)

    output:
    tuple val(species), path("clusters_${params.max_cluster_val}_info.tsv"), emit: max_cluster_info
    tuple val(species), path('clusters_*_info.tsv'),                         emit: cluster_info_file_list
    path("centroid_*_stacked_matrix.tsv")
    path("centroid_${params.max_cluster_val}_lifted_to_${rep_genome}.tsv")

    script:
    """
    cluster_pct_list=\$(echo ${params.centroid_cluster_percents} | tr -d '[],')
    augment_pangenome.py --gene_info_file ${gene_info_tsv} --cluster_thresholds \${cluster_pct_list} --rep_genome ${rep_genome}
    """
}

process CalculateContigLength {
    label 'single_cpu'
    conda "${projectDir}/env/python-bio.yml"
    publishDir "${params.db_path}/pangenomes/${species}/", mode: "copy"

    input:
    tuple val(species), path(fna_list)

    output:
    tuple val(species), path("contigs.len")

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    with open("contigs.len", 'w') as f:
        f.write("genome_id\\tcontig_id\\tcontig_length\\n")
        for genome_fna in "${fna_list}".split(' '):
            gid = genome_fna.replace('.fna','')
            for rec in SeqIO.parse(genome_fna, 'fasta'):
                f.write(f"{gid}\\t{rec.id}\\t{len(rec.seq)}\\n")
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// ANNOTATION RUNNERS
// ─────────────────────────────────────────────────────────────────────────────

process RunEggNog {
    label 'mem_very_high'
    conda "${projectDir}/env/eggnog.yml"
    publishDir "${params.db_path}/pangenomes_annotation/02_eggnog/${species}", mode: "copy"

    input:
    tuple val(species), path(centroid_ffn)
    path eggnog_data

    output:
    tuple val(species), path("${species}.emapper.annotations"), emit: eggnog_output
    path("${species}.emapper*")

    script:
    """
    emapper.py -i ${centroid_ffn} --itype CDS -m diamond --sensmode more-sensitive --data_dir ${eggnog_data} --output ${species} --override --pfam_realign realign --cpu ${task.cpus} --dmnd_db ${eggnog_data}/${params.eggnog_dmnd_db_name}
    """
}

process RunGeNomad {
    label 'mem_high'
    conda "${projectDir}/env/genomad.yml"
    publishDir "${params.db_path}/pangenomes_annotation/01_mge/${species}/${genome}/genomad_output", mode: "copy"

    input:
    tuple val(genome), val(species), path(fna)
    path genomad_db

    output:
    tuple val(genome), val(species), path("${genome}_summary/${genome}_virus_genes.tsv"), path("${genome}_summary/${genome}_plasmid_genes.tsv"), emit: genomad_gene_files
    path("${genome}_*")

    script:
    """
    genomad end-to-end --threads ${task.cpus} --cleanup --enable-score-calibration ${fna} '.' ${genomad_db}
    """
}

process RunResFinder {
    label 'mem_medium'
    conda "${projectDir}/env/resfinder.yml"
    publishDir "${params.db_path}/pangenomes_annotation/01_mge/${species}/${genome}/", mode: "copy"

    input:
    tuple val(genome), val(species), path(fna)
    path res_db
    path dis_db

    output:
    tuple val(genome), val(species), path('resfinder_output/ResFinder_results_tab.txt'), path('mefinder_output/mefinder.csv'), emit: resfinder_output
    path('resfinder_output/*')
    path('mefinder_output/*')

    script:
    """
    mkdir resfinder_output mefinder_output
    cd resfinder_output
    python -m resfinder -ifa ../${fna} -o . -s Other -l ${params.resfinder_min_cov} -t ${params.resfinder_identity_threshold} --acquired -db_res ../${res_db} -d -db_disinf ../${dis_db}
    cd ../mefinder_output
    mefinder find --contig ../${fna} -t ${task.cpus} mefinder
    """
}

process ParsePangenomeAnnotations {
    label 'single_cpu'
    conda "${projectDir}/env/python-bio.yml"
    publishDir "${params.db_path}/pangenomes_annotation/03_processed/${species}/", mode: "copy"

    input:
    tuple val(genome), val(species), path(genes_file), path(resfinder_txt), path(mefinder_csv), path(genomad_virus), path(genomad_plasmid), path(contig_len_file), path(eggnog_results), path(max_cluster_info)

    output:
    tuple val(species), path("genomad_virus/${genome}"), path("genomad_plasmid/${genome}"), path("mefinder/${genome}"), path("resfinder/${genome}"), path("eggnog/${genome}")

    script:
    """
    mkdir genomad_virus genomad_plasmid mefinder resfinder eggnog
    parse_pangenome_annotations.py --species ${species} --genome ${genome} --genes_file ${genes_file} --mefinder_file ${mefinder_csv} --resfinder_file ${resfinder_txt} --genomad_virus_file ${genomad_virus} --genomad_plasmid_file ${genomad_plasmid} --contig_len_file ${contig_len_file} --eggnog_results_file ${eggnog_results} --max_cluster_info_file ${max_cluster_info}
    """
}

process CombinePangenomeAnnotations {
    label 'single_cpu'
    publishDir "${params.db_path}/pangenomes/${species}/annotation", mode: "copy"

    input:
    tuple val(species), path(genomad_virus_files, stageAs: 'genomad_virus/*'), path(genomad_plasmid_files, stageAs: 'genomad_plasmid/*'), path(mefinder_files, stageAs: 'mefinder/*'), path(resfinder_files, stageAs: 'resfinder/*'), path(eggnog_files, stageAs: 'eggnog/*')

    output:
    tuple val(species), path('genomad_virus.tsv'), path('genomad_plasmid.tsv'), path('mefinder.tsv'), path('resfinder.tsv'), path('eggnog.tsv')

    shell:
    '''
    # Header Definitions
    COLS_GENOMAD=('gene_id' 'contig_id' 'start' 'end' 'strand' 'gene_type' 'contig_length' 'start_genomad' 'end_genomad' 'gene' 'annotation_conjscan' 'annotation_amr' 'annotation_accessions' 'annotation_description')
    COLS_MEFINDER=('gene_id' 'contig_id' 'start' 'end' 'strand' 'gene_type' 'contig_length' 'start_mefinder' 'end_mefinder' 'mge_no' 'prediction' 'name' 'type' 'synonyms')
    COLS_RESFINDER=('gene_id' 'contig_id' 'start' 'end' 'strand' 'gene_type' 'contig_length' 'start_resfinder' 'end_resfinder' 'resistance_gene' 'phenotype' 'accession_no')
    COLS_EGGNOG=('#query' 'seed_ortholog' 'evalue' 'score' 'eggNOG_OGs' 'max_annot_lvl' 'COG_category' 'Description' 'Preferred_name' 'GOs' 'EC' 'KEGG_ko' 'KEGG_Pathway' 'KEGG_Module' 'KEGG_Reaction' 'KEGG_rclass' 'BRITE' 'KEGG_TC' 'CAZy' 'BiGG_Reaction' 'PFAMs' 'contig_id' 'start' 'end' 'strand' 'gene_type' 'contig_length')

    # Function to join array with tabs
    function join_by { local d=${1-} f=${2-}; if shift 2; then printf %s "$f" "${@/#/$d}"; fi; }

    # Write headers and append files
    echo -e "$(join_by $'\t' "${COLS_GENOMAD[@]}")" > genomad_virus.tsv
    cat genomad_virus/* >> genomad_virus.tsv

    echo -e "$(join_by $'\t' "${COLS_GENOMAD[@]}")" > genomad_plasmid.tsv
    cat genomad_plasmid/* >> genomad_plasmid.tsv

    echo -e "$(join_by $'\t' "${COLS_MEFINDER[@]}")" > mefinder.tsv
    cat mefinder/* >> mefinder.tsv

    echo -e "$(join_by $'\t' "${COLS_RESFINDER[@]}")" > resfinder.tsv
    cat resfinder/* >> resfinder.tsv

    echo -e "$(join_by $'\t' "${COLS_EGGNOG[@]}")" > eggnog.tsv
    cat eggnog/* >> eggnog.tsv
    '''
}

process EnhancePangenome {
    label 'single_cpu'
    conda "${projectDir}/env/python-bio.yml"
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
    tuple val(species), path(genomad_virus), path(genomad_plasmid), path(mefinder_file), path(resfinder_file), path(eggnog_file), path(genes_info_file), path(cluster_info_file_list, stageAs: 'old/*')

    output:
    tuple val(species), path("genes_annotated.tsv"), path('clusters_*_annot.tsv'), path('clusters_*_info.tsv')

    script:
    """
    cluster_pct_list=\$(echo ${params.centroid_cluster_percents} | tr -d '[],')
    enhance_pangenome.py --cluster_thresholds \${cluster_pct_list} --genomad_virus_file ${genomad_virus} --genomad_plasmid_file ${genomad_plasmid} --mefinder_file ${mefinder_file} --resfinder_file ${resfinder_file} --eggnog_file ${eggnog_file} --genes_info_file ${genes_info_file} --cluster_info_file_list ${cluster_info_file_list}
    """
}

process ComputeChunks {
    label 'single_cpu'
    conda "${projectDir}/env/python-bio.yml"
    publishDir "${params.db_path}/chunks/sites/"

    input:
    tuple val(genome), val(species), path(fna)

    output:
    path("run/chunksize.${params.run_chunk_size}/${species}/${genome}.json")
    path("merge/chunksize.${params.merge_chunk_size}/${species}/${genome}.json")

    script:
    """
    mkdir -p run/chunksize.${params.run_chunk_size}/${species}/ merge/chunksize.${params.merge_chunk_size}/${species}/
    compute_chunks.py --fna ${fna} --genome ${genome} --species ${species} --run_chunk_size ${params.run_chunk_size} --merge_chunk_size ${params.merge_chunk_size}
    """
}