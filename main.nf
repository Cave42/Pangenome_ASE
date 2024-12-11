#!/usr/bin/env nextflow

include { VG_MPMAP } from './modules.nf'
include { RPVG } from './modules.nf'

params.INPUT = false
params.OUTDIR = false


XG = file('bin/1kg_all_af001_gencode100.xg')
GCSA = file('bin/1kg_all_af001_gencode100_index.gcsa')
GCSA_LCP = file('bin/1kg_all_af001_gencode100_index.gcsa.lcp')
DIST = file('bin/1kg_all_af001_gencode100_index.dist')
GBWT = file('bin/1kg_all_af001_gencode100_unidi.gbwt')
GBWT_RI = file('bin/1kg_all_af001_gencode100_unidi.gbwt.ri')
TRANSCRIPT_INFO = file('bin/transcript_info_unidi.txt.gz')


RPVG_SCRIPT = file('bin/rpvg')
RPVG_UNIDI = file('bin/RPVG_UNIDI.sh')
VG_SHELL = file('bin/MpMap.sh')
VG_APP = file('bin/vg_v1.38.0.sif')

/*
XG = file('bin_2/1kg_nonCEU_af001_gencode100.xg')
GCSA = file('bin_2/1kg_nonCEU_af001_gencode100_index.gcsa')
GCSA_LCP = file('bin_2/1kg_nonCEU_af001_gencode100_index.gcsa.lcp')
DIST = file('bin_2/1kg_nonCEU_af001_gencode100_index.dist')
GBWT = file('bin_2/1kg_nonCEU_af001_gencode100.gbwt')
GBWT_RI = file('bin_2/1kg_nonCEU_af001_gencode100.gbwt.ri')
TRANSCRIPT_INFO = file('bin_2/transcript_info_1kg_nonCEU_af001_gencode100_NO_UNIDI.txt.gz')
*/

workflow {
    
    input_read_ch = Channel
        .fromFilePairs("${params.INPUT}*_{1,2}*.gz")
        .ifEmpty { error "Cannot find any FASTQ pairs in ${params.INPUT} ending with .gz" }
        .map { it -> [it[0], it[1][0], it[1][1]]}

    VG_MPMAP (
        input_read_ch,
        XG,
        GCSA,
        DIST,
        GCSA_LCP,
        VG_SHELL,
        VG_APP
    )

    RPVG (
        input_read_ch,
        XG, 
        GBWT,
        TRANSCRIPT_INFO,
        VG_MPMAP.out[0], 
        GBWT_RI,
        RPVG_SCRIPT,
        RPVG_UNIDI
    ) 

}

/*
workflow {

*/

/*
    // Channels for reference files
    xg_file = Channel.fromPath(params.xg_file)
    gcsa_file = Channel.fromPath(params.gcsa_file)
    dist_file = Channel.fromPath(params.dist_file)
    gbwt_file = Channel.fromPath(params.gbwt_file)
    transcript_info = Channel.fromPath(params.transcript_info)
    
    // Step 1: VG_MPMAP
    gamp_files = fastq_pairs
                    .map { pair -> tuple(pair[0].baseName.replaceAll(/_1\.fastq\.gz$/, ""), pair[0], pair[1]) }
                    .set { fastq_pairs }
                    | VG_MPMAP(xg_file, gcsa_file, dist_file)
    
    // Step 2: RPVG, using output of VG_MPMAP as input
    gamp_files | RPVG(xg_file, gbwt_file, transcript_info)
    */