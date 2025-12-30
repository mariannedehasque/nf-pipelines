nextflow.enable.dsl=2

// --- Default Parameters ---
params.samplesheet = "${projectDir}/inputfiles/samplesheet.csv"
params.contigs     = "${projectDir}/inputfiles/contigs.txt"
params.outdir       = "${projectDir}/results" 
params.reference    = "${projectDir}/data/reference/reference.denovoSSL.Tzo20k.fasta"
params.bed_file     = "${projectDir}/data/reference/reference.denovoSSL.Tzo20k.repma.angsd.txt"
params.species      = "Tzo"
params.maxdepth    = 10000
params.minind      = 3

// import modules
include { ANGSD_GL_ALL; ANGSD_COLLECT_OUTPUT; ANGSD_EXTRACT_SITES } from './modules/angsd_gl'
include { COLLECT_BAM_ALL; COLLECT_BAM_POP } from './modules/collect_bam'
include { PCANGSD; PLOT_PCANGSD } from './modules/pcangsd'


// --- Input Channels ---

// Read in samplesheet and collect metadata

ch_samples = channel
    .fromPath(params.samplesheet, checkIfExists: true)
    .splitCsv(header: true)
    .map { row -> [row.sample,row.pop,row.bam]}

ch_pop_bams = ch_samples
    .groupTuple(by: 1)
    .map { samples, pop, bams -> tuple(pop, bams)}

ch_all_bams = ch_samples
    .map { sample, pop, bam -> bam }
    .collect()

ch_contigs = Channel
    .fromPath("${params.contigs}", checkIfExists: true)
    .splitText()
    .map { line -> line.trim() }


// --- Workflow ---

workflow {
    // Create bamlist files
    ch_samplesheet_file = Channel.fromPath(params.samplesheet)
    ch_bamlist = COLLECT_BAM_ALL(ch_all_bams)
    //ch_bamlist_pop = COLLECT_BAM_POP(ch_pop_bams)
    
    // ANGSD genotype likelihoods
    ch_genotypes = ANGSD_GL_ALL(ch_contigs, ch_bamlist)
    ch_mafs = ch_genotypes.mafs.map{contig, file -> file }.collect()
    ch_beagle = ch_genotypes.beagle.map{contig, file -> file }.collect()
    ch_collected = ANGSD_COLLECT_OUTPUT(ch_mafs, ch_beagle)
    ch_sites = ANGSD_EXTRACT_SITES(ch_collected.all_mafs)
    
    // PCAngsd
    //ch_pcangsd_cov = PCANGSD(ch_collected.all_beagle)
    //PLOT_PCANGSD(ch_pcangsd_cov, ch_samplesheet_file)
}