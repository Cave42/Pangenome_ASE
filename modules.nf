process VG_MPMAP {

    //errorStrategy 'ignore'
    //memory { task.attempt <= 1 ? '72 GB' : '300 GB' }
    //maxRetries 2
    //errorStrategy { task.attempt > 2 ? 'ignore' : 'retry' }
    //errorStrategy 'ignore'


    publishDir params.OUTDIR
    
    input:
    tuple val(base), file(R1), file(R2)
    file XG
    file GCSA
    file DIST
    file GCSA_LCP
    file VG_SHELL
    file VG_APP

    output:
    file("mpmap_1kgALL_${base}.gamp") 

    container 'bin/vg_v1.38.0.sif'
    
    script:

    """

    vg mpmap -n rna -t ${task.cpus} -x ${XG} -g ${GCSA} -d ${DIST} -f ${R1} -f ${R2} > mpmap_1kgALL_${base}.gamp

    """
}

process RPVG {

    //maxForks 5

    errorStrategy 'ignore'

    publishDir params.OUTDIR

    input:
    tuple val(base), file(R1), file(R2)
    file XG
    file GBWT
    file TRANSCRIPT_INFO 
    file GAMP
    file GBWT_RI
    file RPVG_SCRIPT
    file RPVG_UNIDI

    output:

    file('*results_joint.txt')
    file('*results.txt')

    shell:

    """

    srun --account=cpu-gdml-sponsored --partition=cpu-core-sponsored --cpus-per-task 32 --mem=124G --time=24:00:00 ./RPVG_UNIDI.sh ${XG} \
                  ${GBWT} \
                  ${GAMP} \
                  ${base}_rpvg_results \
                  ${TRANSCRIPT_INFO}

    """
}

/*
    # srun --account=cpu-gdml-sponsored --partition=cpu-core-sponsored --cpus-per-task 32 --mem=124G --time=24:00:00 --mail-user cave42@uw.edu ./MpMap.sh ${XG} \
    #              ${GCSA} \
    #              ${DIST} \
    #              ${R1} \
    #              ${R2} \
    #              mpmap_1kgALL_${base}.gamp
*/

//vg mpmap -n rna -t ${task.cpus} -x ${XG} -g ${GCSA} -d ${DIST} -f ${R1} -f ${R2} > mpmap_1kgALL_${base}.gamp
//rpvg -t ${task.cpus} -g ${XG} -p ${GBWT} -a ${GAMP} -o ${base}_rpvg_results -i haplotype-transcripts -f ${TRANSCRIPT_INFO}
//srun --account=cpu-gdml-sponsored --partition=cpu-core-sponsored --cpus-per-task 32 --mem=124G --time=24:00:00 ./RPVG.sh