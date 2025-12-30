nextflow.enable.dsl=2

// --- Default Parameters ---
params.samplesheet = "${projectDir}/inputfiles/samplesheet.csv"
params.contigs     = "${projectDir}/inputfiles/contigs.txt"
params.outdir       = "${projectDir}/results" 
params.reference    = "${projectDir}/data/reference/reference.denovoSSL.Tzo20k.fasta"
params.bed_file     = "${projectDir}/data/reference/reference.denovoSSL.Tzo20k.repma.angsd.txt"
params.species      = "Tzo"


// --- Input Channel ---

// Read in samplesheet and collect metadata
ch_samples = channel
    .fromPath(params.samplesheet, checkIfExists: true)
    .splitCsv(header: true)
    .map { row -> [row.sample,row.pop,row.bam]}

ch_pop_bams = ch_samples
    .groupTuple(by: 1)
    .map { samples, pop, bams -> tuple(pop, bams)}
    .view()

ch_all_bams = ch_samples
    .map { sample, pop, bam -> bam }
    .collect()

ch_contigs = Channel
    .fromPath("${params.contigs}", checkIfExists: true)
    .splitText()
    .map { line -> line.trim() }