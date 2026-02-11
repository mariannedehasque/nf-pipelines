## Overview

This pipeline trims fastq files and maps using ancient DNA-specific settings. It is written in Nextflow DSL2 and designed to be run on Old Dominion WAHAB cluster. The pipeline uses the same symbolic filename structure as the GenErode pipeline.

## Folder structure

## Parameters

### Inputfiles

The pipeline requires the following inputfiles:

* Reference (fasta or fna format)
* Reference index files (.bwt, .ann, .sa, .pac, .ann, .amb )
* BED file with masked repeat regions
* File with the name of all fastq files to be processed 
* Directory containing all fastq files

All inputfiles can be copied or generated from the GenErode directory.

```bash
# Activate bash shell

bash

# Change directory to the nextflow pipeline
cd ./nf-trim-merged-unmerged

# Create new directories
mkdir data
mkdir ./data/reference
mkdir ./data/symlinks
mkdir inputfiles

# Create softlinks to the raw fastq files
ln -s /Generode/data/raw_reads_symlinks/modern/*fastq.gz ./data/symlinks

# Create fastq filenames file. Manually adjust the file if necessary (e.g. if not all samples from GenErode are to be used)

ls ./data/symlinks/*fastq.gz | xargs -n1 basename | cut -d "_" -f1,2,3 | uniq > ./inputfiles/fastq_filenames.txt

# Create softlinks to reference and repma bed file
# Adjust the path to the reference if necessary

ln -s /Generode/reference/<reference>.fasta ./data/reference/
ln -s /Generode/reference/<reference>.fasta.* ./data/reference/
ln -s /Generode/reference/<reference>.repma.bed ./data/reference/

```

### Configuration

In the `main.nf` file, adjust the parameters. If you used the same structure as outlined above, only the name of the reference should be added.

Params.trimlength will depend on the average read length of historical samples. I advice to keep the trim length < 89

```bash
// --- Default Parameters ---
params.samples_file = "${projectDir}/inputfiles/fastq_filenames.txt"
params.indir        = "${projectDir}/data/symlinks"
params.outdir       = "${projectDir}/results" 
params.reference    = "${projectDir}/data/reference/<reference>.fasta"
params.bed_file     = "${projectDir}/data/reference/<reference>.repma.bed"
params.split_script = "${projectDir}/scripts/split_reads.sh"
params.rmdup_script = "${projectDir}/scripts/samremovedup.py"
params.amber_script = "/home/mdehasqu/TOOLS/AMBER/AMBER" // Ignore this.
params.bwa_threads  = 4
params.bam_q        = 1 // Mapping quality. Currently set to 1 simply to remove unmapped reads. 
params.trimlength   = 85
```

Since the `AMBER` tool requires special software and a unique conda environment on WAHAB, we will just skip this step for now.

At the bottom of the `main.nf` file, hash out the AMBER rules.

```bash
    // Run AMBER
    // AMBER_PREP(INDEX_REALIGNED.out)
    // AMBER(AMBER_PREP.out)
```

## To Run

I like to run nextflow pipelines in tmux so that I can keep it running in the background even when logging of.

To create a new tmux screen:

```bash
tmux new -s nextflow
```

!Good to know! Leaving a tmux sessions is notoriously difficult. To leave and reopen, type the following commands:

```bash
# To leave, literraly press these keys
Ctrl + B, followed by D

# To re-enter a session
tmux a -t nextflow

```

Then start the nextflow pipeline as follows (while in tmux):

```bash
# Load the nextflow module
module load container_env
module load nextflow

# Start the run
nextflow run main.nf -profile standard -resume

```

If everything went well, the pipeline will run and submit jobs to the queue. 
On the first run, a new conda environment will be created. This can take some time.