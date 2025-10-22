#!/usr/bin/env nextflow

include { mapReads } from './modules/local/mapReads.nf'
include { samIndex } from './modules/local/samIndex.nf'
include { bamCoverage } from './modules/local/bamCoverage.nf'
include { readsCount } from './modules/local/readsCount.nf'
//include { collectReadCounts } from './modules/local/collectReadCounts.nf'

params.inputFile='inputFile_main_placeholder'

workflow{
    genomePath_ch=Channel.fromPath(params.inputFile)
                         .splitCsv(header: true)
                         .map { item -> item.Genome }

    readsPath_ch=Channel.fromPath(params.inputFile)
                        .splitCsv(header: true)
                        .map { item -> item.Reads }


    genomeName_ch=Channel.fromPath(params.inputFile)
                         .splitCsv(header:true)
                         .map { item -> item.Genome_Name }

    distinctReads_ch=Channel.fromPath(params.inputFile)
                            .splitCsv(header: true)
                            .map { item -> item.Reads }
                            .flatten()
                            .distinct()

    readsCount(distinctReads_ch)

    mappedBam_ch = mapReads(genomePath_ch,readsPath_ch,genomeName_ch)

    samIndex(mappedBam_ch)

    bamCoverage(mappedBam_ch)
}