process ANGSD_DIVERSITY {
    tag { pop }
    publishDir "${params.outdir}/angsd_pop_theta", mode: 'copy'

    input:
    tuple val(pop), path(saf), path(saf_idx), path(pos)

    output:
    path "${pop}.pestPG", emit: pestPG
    path "${pop}.thetas.idx", emit: thetas_idx
    path "${pop}.thetas.gz", emit: thetas
    path "${pop}.sfs", emit: sfs

    script:
    """
        realSFS ${pop}.saf.idx -P 8 -fold 1 > ${pop}.sfs
        realSFS saf2theta ${pop}.saf.idx -sfs ${pop}.sfs -fold 1 -P 8 -outname ${pop}
        thetaStat do_stat ${pop}.thetas.idx -outnames ${pop}
    """
}