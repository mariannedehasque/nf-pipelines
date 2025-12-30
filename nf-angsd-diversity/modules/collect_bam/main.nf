process COLLECT_BAM_ALL {

    publishDir "${params.outdir}/inputfiles", mode: 'copy'

    input:
    val bams

    output:
    path "bamlist.txt"

    script:
    """
    printf '%s\\n' ${bams.join(' ')} | head -c -1> bamlist.txt
    """
}

process COLLECT_BAM_POP {

    tag { pop }
    publishDir "${params.outdir}/inputfiles", mode: 'copy'

    input:
    tuple val(pop), val(bams)

    output:
    path "${pop}.bamlist.txt"

    script:
    """
    printf '%s\\n' ${bams.join(' ')} | head -c -1 > ${pop}.bamlist.txt
    """
}