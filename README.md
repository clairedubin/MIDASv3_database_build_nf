## Overview

This Nextflow pipeline is designed to create a custom database for use with the MIDAS profiling tool. 

## Installation

1. [Install nextflow](https://www.nextflow.io/docs/latest/install.html). This pipeline has been tested on Nextflow versions 23.10.0 and 24.10.5.
2. [Update nextflow.config](https://www.nextflow.io/docs/latest/config.html) for your system or HPC. The included nextflow.config file is specific to the developer's HPC (SGE) and will not work for other systems.
3. Install the Eggnog database.
4. Test nextflow installation.

## Inputs 

Genome assemblies must already be clustered into species or species-level genome bins.

A TSV with five columns:
- `genome`: the name of the genome
- `species`: 6 digit identifier for the species - these can be randomly assigned.
- `representative`: the representative genome from the species cluster
- `genome_is_representative`: whether the genome is the representative genome for the species cluster
- `fasta_path`: full path to the genome assembly file

Example genomes.tsv:

| genome            | species  | representative        | genome_is_representative | fasta_path                                    |
|-------------------|----------|-----------------------|--------------------------|-------------------------------------------------|
| GCA_900552055.1   | 117086   | GCA_900552055.1       | 1                       | /Users/myname/MIDAS3_nextflow/testing/inputs/genomes/GCA_900552055.1.fasta |
| GCF_900752885.1   | 117086   | GCA_900552055.1       | 0                 | /Users/myname/MIDAS3_nextflow/testing/inputs/genomes/GCF_900752885.1.fasta |
| GCA_007120565.1   | 117088   | GCA_007120565.1       | 1                 | /Users/myname/MIDAS3_nextflow/testing/inputs/genomes/GCA_007120565.1.fasta |
| GCA_007135585.1   | 117088   | GCA_007120565.1       | 0                 | /Users/myname/MIDAS3_nextflow/testing/inputs/genomes/GCA_007135585.1.fasta |

## Required parameters

The following parameters must be included in your Nextflow run command or via a `params.json` file:

| Parameter                   | Description                          | Example Value                                        |
| --------------------------- | ------------------------------------ | ----------------------------------------------------- |
| `--genomes_tsv_path`        | Path to genomes file    | `genomes.tsv`   |
| `eggnog_db_dir`             | Path to Eggnog database directory    | `/wynton/group/sirota/clairedubin/databases/eggnog`   |
| `eggnog_dmnd_db_name`       | Eggnog Diamond database name         | `eggnog_proteins.dmnd`                                |
| `eggnog_conda_dir`          | Path to Eggnog Conda environment     | `/wynton/protected/home/sirota/clairedubin/anaconda3/envs/eggnog` |
| `genomad_db_dir`            | Path to Genomad database directory   | `/wynton/protected/home/sirota/clairedubin/databases/genomad_db_v1.5/genomad_db` |
| `genomad_conda_dir`         | Path to Genomad Conda environment    | `/wynton/protected/home/sirota/clairedubin/anaconda3/envs/genomad` |
| `resfinder_env_dir`          | Path to ResFinder environment        | `/wynton/protected/home/sirota/clairedubin/envs/resfinder_env` |
| `resfinder_db_dir`           | Path to ResFinder databases          | `/wynton/protected/home/sirota/clairedubin/databases/resfinder_dbs` |
| `blastn_dir`                | Path to BLASTN executable            | `/wynton/protected/home/sirota/clairedubin/bin/ncbi-blast-2.14.1+/bin` |
| `git_dir`                   | Path to Git executable               | `/wynton/protected/home/sirota/clairedubin/bin/git-2.39.5`           |


## Optional Parameters

The following parameters can be configured in your Nextflow run command or in the `params.json` file:

| Parameter                   | Description                          | Default Value                                         |
| --------------------------- | ------------------------------------ | ----------------------------------------------------- |
| `--db_name`                   | Name of the database                 | `nextflow_db`                                         |
| `--db_output_dir`             | Directory for database output        | `.` (current directory)                               |
| `--centroid_cluster_percents` | Cluster percentages for centroids    | `[99, 95, 90, 85, 80, 75]`                           |
| `--run_chunk_size`            | Chunk size for running processes     | `1000000`                                             |
| `--merge_chunk_size`          | Chunk size for merging processes     | `500000`                                              |
| `--resfinder_min_cov`         | Minimum coverage for ResFinder       | `0.6`                                                 |
| `--resfinder_identity_threshold` | Identity threshold for ResFinder   | `0.8`                                                 |
| `--hmmsearch_min_cov`         | Minimum coverage for HMM search      | `0.0`                                                 |
| `--hmmsearch_max_evalue`      | Maximum e-value for HMM search       | `1e-5`                                                |
| `--marker_set`                | Marker set for HMM models            | `phyeco`                                              |


The following nextflow command line flags can also be used. Additional flags and information can be found in the [Nextflow documentation](https://www.nextflow.io/docs/stable/cli.html).

| Parameter                   | Description                          | Default Value                                         |
| --------------------------- | ------------------------------------ | ----------------------------------------------------- |
| `-w`                | Directory for scratch data            | `./work`                                                |
| `-c`                | Configuration file            | `nextflow.config`                                              |
| `-params-file`                | Parameters file            | `params.json`                                              |
| `-resume`                | Resume previously running job            |                                |
| `-h` or   `--help`        | Show help message and exit            |                                |



## Usage

To run the pipeline:

```bash
nextflow run main.nf --db_output_dir /path/to/db/output/ --db_name example_db
```

You can also use a configuration file:

```bash
nextflow run main.nf -c nextflow.config
```


## Example data

Example data is available in the testing folder. To build a database from this data, run:
```bash
cd ~/MIDAS3_nextflow #location of this nextflow pipeline
nextflow run main.nf -c nextflow.config -params-file params.json --genomes_tsv_path testing/inputs/genomes.tsv --db_output_dir testing -w testing/work```
```

## Citation

If you use this tool, please cite:
TODO: MIDAS paper
TODO: My paper!

