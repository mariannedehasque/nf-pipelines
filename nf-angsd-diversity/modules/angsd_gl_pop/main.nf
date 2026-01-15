process ANGSD_GL_POP {
    tag { pop }
    publishDir "${params.outdir}/angsd_pop", mode: 'copy'

    input:
    tuple val(pop), path(bamlist)
    path snps
    path bin
    path idx
    path regions

    output:
    tuple val(pop), path("${pop}.saf.gz"), path("${pop}.saf.idx"), path("${pop}.saf.pos.gz"), emit: saf_files
    path "${pop}.beagle.gz", emit: beagle
    path "${pop}.mafs.gz", emit: mafs

    script:
    """
        angsd -bam ${bamlist} \
        -doSaf 1 \
        -GL 1 \
        -doGlf 2 \
        -doMajorMinor 4 \
        -doMaf 1 \
        -minMapQ 25 -minQ 30 \
        -doCounts 1 -doDepth 1 -dumpCounts 1 \
        -uniqueOnly 1 -remove_bads 1 \
        -P ${task.cpus} \
        -ref ${params.reference} \
        -anc ${params.reference} \
        -sites ${snps} \
        -rf ${regions} \
        -out ${pop}
    """
}