# MIDAS3 Database Build Pipeline

## Overview

This Nextflow pipeline is designed to create a custom database from a set of prokaryote genome assemblies for use with MIDAS (version 3). 

---

## Installation

### 1. Requirements

- **Nextflow:** Install [Nextflow](https://www.nextflow.io/). This pipeline has been tested on versions 23.10.0 and 24.10.5.
- **Conda:** Install [Conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda or Anaconda) for environment management.

### 2. Prebuild Conda Environments

To save time during execution and ensure stability on HPC clusters, run the prebuild script. This creates all necessary environments in a `.conda_envs/` directory and generates a configuration file so Nextflow knows where to find them.
```bash
bash bin/prebuild_conda_envs.sh
```

> **Note:** This modifies `conf/conda_envs.config`. If you want to revert to on-the-fly environment creation from YML files, re-download that file from Github.

### 3. HPC Configuration

The included `nextflow.config` is pre-configured for SGE on UCSF's Wynton cluster. Update the `process` and `executor` settings in `nextflow.config` to match your specific system or HPC (e.g., Slurm, LSF).

---

## Database Setup

This pipeline handles the installation of the EggNOG, geNomad, and ResFinder databases automatically.

- Specify the desired installation paths for these databases in your `params.json`.
- If the specified directories do not exist, the pipeline will automatically trigger download processes (run locally).
  
---

## Inputs

Genome assemblies must already be clustered into species or species-level genome bins. Provide a TSV file with the following five columns:

| Column | Description |
|--------|-------------|
| `genome` | The name of the genome. |
| `species` | A unique identifier for the species cluster (e.g., a 6-digit ID). |
| `representative` | The name of the representative genome for that cluster. |
| `genome_is_representative` | `1` if the genome is the representative, `0` otherwise. |
| `fasta_path` | The full path to the genome assembly file. |

**Example `genomes.tsv`:**

| genome | species | representative | genome_is_representative | fasta_path |
|--------|---------|----------------|--------------------------|------------|
| GCA_900552055.1 | 117086 | GCA_900552055.1 | 1 | /path/to/genomes/GCA_900552055.1.fasta |
| GCF_900752885.1 | 117086 | GCA_900552055.1 | 0 | /path/to/genomes/GCF_900752885.1.fasta |

---

## Parameters

### Required Parameters (`params.json`)

| Parameter | Description |
|-----------|-------------|
| `eggnog_db_dir` | Directory where EggNOG data is stored/downloaded. |
| `genomad_db_dir` | Directory where geNomad data is stored/downloaded. |
| `resfinder_db_dir` | Directory where ResFinder/PointFinder data is stored/downloaded. |
| `eggnog_dmnd_db_name` | Name of the EggNOG Diamond database (usually `eggnog_proteins.dmnd`). |
| `marker_set` | Marker set for HMM models (default: `phyeco`). |

### Optional Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--db_name` | Name of the output database. | `nextflow_db` |
| `--db_output_dir` | Parent directory for the database output. | `.` |
| `--centroid_cluster_percents` | Identity thresholds for clustering. MIDAS expects databases to have the default thresholds. | `[99, 95, 90, 85, 80, 75]` |

---

## Usage

To run the pipeline with a parameters file:
```bash
nextflow run main.nf -c nextflow.config -params-file params.json --genomes_tsv_path genomes.tsv
```

### Resume Execution

If a run fails or is interrupted, use the `-resume` flag to pick up where it left off:
```bash
nextflow run main.nf -c nextflow.config -params-file params.json -resume
```

### Example Test Run
```bash

bash bin/prebuild_conda_envs.sh

nextflow run main.nf \
  -c nextflow.config \
  -params-file params.json \
  --genomes_tsv_path testing/inputs/genomes.tsv \
  --db_output_dir testing \
  -w testing/work
```

---

## Citation

If you use this tool, please cite:

- Dubin CA, Zhao C, Pollard KS, Oskotsky T, Golob JL, Sirota M. Expanding vaginal microbiome pangenomes via a custom MIDAS database reveals *Lactobacillus crispatus* accessory genes associated with cervical dysplasia. *bioRxiv*. 2025. doi:[10.1101/2025.09.11.675634](https://doi.org/10.1101/2025.09.11.675634)

- Smith BJ, Zhao C, Dubinkina V, Jin X, Zahavi L, Shoer S, Moltzau-Anderson J, Segal E, Pollard KS. Accurate estimation of intraspecific microbial gene content variation in metagenomic data with MIDAS v3 and StrainPGC. *Genome Res*. 2025. doi:[10.1101/gr.279543.124](https://doi.org/10.1101/gr.279543.124)
