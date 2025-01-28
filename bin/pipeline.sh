# Adapted from Chunyu Zhao, originally from 2023-06-14
#!/usr/bin/bash

#Arguments
species_id=$1
genes_info=$2
centroids_clean_ffn=$3
centroids_ambig_ffn=$4
genes_ffn=$5
genes_len=$6
total_threads=$7
total_mem=$8
script_dir=$9

# Directories and Variables
base_dir=$PWD
cluster_threshold=$(echo "$centroids_clean_ffn" | grep -oP '(?<=centroids\.)\d+(?=\.clean\.ffn)')

out_dir="$PWD/temp/cdhit"
members_dir="${out_dir}/step1_members"
vsearch_dir="${out_dir}/step2_vsearch"
cdhit_dir="${out_dir}/step3_cdhit"
info_dir="${out_dir}/step4_info"

mkdir -p "$out_dir" "$members_dir" "$vsearch_dir" "$cdhit_dir" "$info_dir"

# File paths
gene_info_vsearch="${info_dir}/gene_info_vsearch.tsv"
vsearch_centroids_ffn="${info_dir}/vsearch_centroids.ffn"
cdhit_centroids_ffn="${info_dir}/cdhit_centroids.ffn"
cdhit_centroids_tsv="${info_dir}/cdhit_centroids.tsv"
gene_info_cdhit="${info_dir}/gene_info_cdhit.tsv"

# Step 1: Handle Ambiguous Centroids
if [[ -s "$centroids_ambig_ffn" ]]; then
  if [[ ! -e "$gene_info_vsearch" ]]; then
    ambiguous_list="${members_dir}/list_of_ambiguous_centroids"
    grep ">" "$centroids_ambig_ffn" | awk '{sub(/>/, ""); print}' > "$ambiguous_list"

    cat "$ambiguous_list" | \
      xargs -Ixx -P "$total_threads" bash -c "bash $script_dir/get_members.sh xx $genes_ffn $genes_info $members_dir/xx.mems.ffn"

    multimember_centroids="${members_dir}/centroids_with_multimembers"
    find "$members_dir" -name '*mems.ffn' -size +0c | awk '{sub(/\.mems\.ffn/, ""); print $NF}' > "$multimember_centroids"

    cat "$multimember_centroids" | \
      xargs -Ixx -P "$total_threads" bash -c "vsearch --cluster_fast $members_dir/xx.mems.ffn --threads 1 --quiet --id 0.99 --centroids $vsearch_dir/xx.centroids -uc $vsearch_dir/xx.clusters"

    cat "$multimember_centroids" | \
      xargs -Ixx -P "$total_threads" bash -c "bash $script_dir/gather_geneinfo.sh $vsearch_dir/xx.clusters $vsearch_dir/xx.geneinfo"

    gene_info_add="${info_dir}/vsearch_gene_info_add.tsv"
    cat "$vsearch_dir"/*.geneinfo | tr ' ' '\t' > "$gene_info_add"
    cut -f2 "$gene_info_add" | sort | uniq > "$vsearch_dir/list_of_vsearch_centroids"
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$0]; next} !($2 in a)' "$ambiguous_list" "$genes_info" > "$vsearch_dir/vsearch_gene_info_keep.tsv"

    cat "$vsearch_dir/vsearch_gene_info_keep.tsv" "$gene_info_add" > "$gene_info_vsearch"
  fi

  if [[ ! -e "$vsearch_centroids_ffn" ]]; then
    seqkit grep -w 0 -f "$vsearch_dir/list_of_vsearch_centroids" "$genes_ffn" > "$vsearch_dir/centroids_add.ffn"
    cat "$centroids_clean_ffn" "$vsearch_dir/centroids_add.ffn" > "$vsearch_centroids_ffn"
  fi
else
  echo "No ambiguous centroids_${cluster_threshold}"
  cp "$genes_info" "$gene_info_vsearch"
  cp "$centroids_clean_ffn" "$vsearch_centroids_ffn"
fi

# Step 2: CD-HIT Processing
if [[ ! -e "$cdhit_centroids_tsv" ]]; then
  pushd "$cdhit_dir"
  cp "$vsearch_centroids_ffn" "vsearch_centroids.ffn"

  awk '/^>/{$0=">g"++i; }1' vsearch_centroids.ffn > vsearch_centroids_renamed.ffn
  paste <(grep "^>" vsearch_centroids.ffn | cut -c 2-) <(grep "^>" vsearch_centroids_renamed.ffn | cut -c 2-) > mapping.txt

  cd-hit-est -i vsearch_centroids_renamed.ffn -c 1 -T "$total_threads" -aS 0.9 -G 0 -g 1 -AS 180 -M "$total_mem" -o cdhit_centroids.ffn

  awk '{ if ($0 ~ /^>/) {print $0} else {print toupper($0)}}' cdhit_centroids.ffn > cdhit_centroids_upper.ffn
  bash "$script_dir/parse_cdhit_cluster.sh" cdhit_centroids.ffn.clstr cdhit_centroids.ffn.clstr.tsv

  awk -v OFS='\t' '$2 == 0 {print $1, $4, $4}' cdhit_centroids.ffn.clstr.tsv > new_centroids.tsv
  awk -v OFS='\t' '$2 != 0 {print $1, $4}' cdhit_centroids.ffn.clstr.tsv > new_members.tsv
  join -t $'\t' <(sort -k1 new_members.tsv) <(sort -k1 new_centroids.tsv | cut -f1-2) > new_members_w_centroids.tsv

  cat new_centroids.tsv new_members_w_centroids.tsv | sort -k1,2n > cdhit_centroids.tsv

  awk 'BEGIN{while(getline <"mapping.txt") dict[">"$2] = ">"$1} {if ($0 in dict) print dict[$0]; else print $0}' cdhit_centroids_upper.ffn > "$cdhit_centroids_ffn"
  awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$2]=$1; next} $2 in a {$2=a[$2]} 1' mapping.txt cdhit_centroids.tsv > cdhit_centroids_c2.tsv
  awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$2]=$1; next} $3 in a {$3=a[$3]} 1' mapping.txt cdhit_centroids_c2.tsv > "$cdhit_centroids_tsv"

  popd 
fi

# Step 3: Gene Information Mapping
if [[ ! -e "$gene_info_cdhit" ]]; then
  join -t $'\t' <(awk 'BEGIN{OFS=FS="\t"}{print $2, $1}' "$gene_info_vsearch" | sort -k1,1) <(awk 'BEGIN{OFS=FS="\t"}{print $2, $3}' "$cdhit_centroids_tsv" | sort -k1,1) > "${info_dir}/cdhit_gene_centroid_mapping.tsv"
  awk -v OFS='\t' '{print $2, $3}' "${info_dir}/cdhit_gene_centroid_mapping.tsv" > "$gene_info_cdhit"
  cut -f1 "$gene_info_cdhit" > "${info_dir}/list_of_cdhit_genes"
fi

# Step 4: Final Outputs
awk '{ if ($0 ~ /^>/) {print $0} else {print toupper($0)}}' "$cdhit_centroids_ffn" > "${out_dir}/centroids.${cluster_threshold}.ffn"
seqkit grep -w 0 -f "${info_dir}/list_of_cdhit_genes" "$genes_ffn" > "${out_dir}/genes.ffn"
cp "$gene_info_cdhit" "${out_dir}/gene_info.txt"
grep -Fwf <(cut -f1 "$gene_info_cdhit") "$genes_len" > "${out_dir}/genes.len"
