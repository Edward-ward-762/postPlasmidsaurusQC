# postPlasmisaurusQC
Nextflow workflow used for processing fastq files from Plasmidsaurus's premium PCR offering, but the general functionality will apply to any fastq file.

The workflow currently calculates the number of reads in the the input fastq files, and summarises them in a single csv, and will aligns the fastq files against a genome of choice. From the genome aligned reads a bedgraph is file is calculated, showing where in the genome your reads map, allowing you to check for undesired mis-priming activity.

## Installation guide
The workflow should be entirely self-contained apart from two dependencies:
 - Nextflow
 - docker

The workflow will work using docker containers, I haven't tested other container options, but I wouldn't expect them to work without a local installation of samtools.

## Usage guide

### Input Files

#### Input csv

The inputFile.csv has 3 required columns and is laid out as follows:

|sample_id |fastq_path           |genome_path        |
|----------|---------------------|-------------------|
|sample_1  |/path/to/fastq.fastq |/path/to/genome.fa |

- **sample_id:**
    - **Description:** A unique name for each sample (row)

- **fastq_path:**
    - **Description:** A file path to your fastq to analyse

- **genome_path:**
    - **Description:** A file path to the genome you'd like to map against. This can be left blank if you don't want to map against a genome. ie: plasmid sequencing

#### Bash script

Included in the git repository is a minimal example bash script to run the pipeline. It has two functions. It will pull the latest version of the pipeline from the main repository. It will then run the pipeline with the docker profile as standard, with -resume, so it will start from where it stops if something is wrong (file path not correct). The inputFile parameter will need to be changed, to point at your own input file, detailed above.

Parameters you can change are outlined below.

### Workflow parameters

#### Output Options

- **`outdir`**  
  _Default:_ `./results`  
  **Description:** Directory where the pipeline's results will be saved.

- **`tracedir`**  
  _Default:_ `.${params.outdir}/pipeline_info`  
  **Description:** Directory where the pipeline info (run reports, software versions) will be saved.

- **`publish_dir_mode`**  
  _Default:_ `copy`  
  **Description:** Mode Nextflow will use when handling output of processes.

  **Known possible values:**
  - `copy` - default value
  - `copyNoFollow`
  - `link`
  - `move`
  - `rellink`
  - `symlink`
  For more information regarding how each option will function please consult: https://www.nextflow.io/docs/latest/reference/process.html#publishdir

- **`debug`**  
  _Default:_ `false`  
  **Description:** Enables debug mode if set to true.

#### Max Resource Options

- **`max_memory`**  
  _Default:_ `128.GB`  
  **Description:** Maximum memory allocation for the pipeline **per process**. Processes may request less.

- **`max_cpus`**  
  _Default:_ `16`  
  **Description:** Maximum number of CPUs that can be used **per process**. Processes may request less.

- **`max_time`**  
  _Default:_ `240.h`  
  **Description:** Maximum time allocation for the pipeline **per process**. Processes may request less.

#### Workflow parameters

- **`inputFile`**
  _Default:_ `inputFile_main_placeholder`
  **Description:** File path to input csv file

- **`min_bam_size`**
  _Default:_ `500`  
  **Description:** Minimum bam file size in bytes for an empty bam file. Used in filtering out empty bam files.

- **`bedgraph_bin_size`**
  _Default:_ `10000`
  **Description:** Bin size used when creating bedgraph file. Genome will be split into bins of this size, and the number of reads counted in each bin.
