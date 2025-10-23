#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    post plasmidsaurus QC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { mapReads } from './modules/local/mapReads.nf'
include { samIndex } from './modules/local/samIndex.nf'
include { bamCoverage } from './modules/local/bamCoverage.nf'
include { readsCount } from './modules/local/readsCount.nf'
include { collectReadCounts } from './modules/local/collectReadCounts.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow{

    //
    // ****************************
    //
    // SECTION: Find number of reads in distinct fastqs
    //
    // ****************************
    //

    //
    // CHANNEL: Find distinct fastq file paths
    //

    ch_distinctReads = Channel.fromPath(params.inputFile)
        .splitCsv(header: true)
        .map { item -> item.Reads }
        .flatten()
        .distinct()

    //
    // MODULE: Find number of reads in fastq
    //

    readsCount(
        ch_distinctReads
    )

    //
    // MODULE: Pool distinct reads into a single csv file
    //

    collectReadCounts(
        readsCount.out
            .collect()
    )

    //
    // ****************************
    //
    // SECTION: Map all read in fqs to genome
    //
    // ****************************
    //

    //
    // CHANNEL: create channel to replace those below
    //

    ch_data = Channel.fromPath(params.inputFile)
                .splitCsv(header: true)
                .map { row ->
                    [[id: row.sample_id, genome_path: row.genome_path], row.fastq_path]
                }

    //
    // CHANNEL: create versions channel
    //

    ch_versions = Channel.empty()

    //
    // CHANNEL: create genomePath channel from genome file path
    //

    genomePath_ch=Channel.fromPath(params.inputFile)
                         .splitCsv(header: true)
                         .map { item -> item.Genome }

    //
    // CHANNEL: create readsPath channel from fq file path
    //

    readsPath_ch=Channel.fromPath(params.inputFile)
                        .splitCsv(header: true)
                        .map { item -> item.Reads }

    //
    // CHANNEL: create genomeName channel from input custom genome name
    //

    genomeName_ch=Channel.fromPath(params.inputFile)
                         .splitCsv(header:true)
                         .map { item -> item.Genome_Name }

    //
    // MODULE: map fq reads to genome
    //

    mappedBam_ch = mapReads(genomePath_ch,readsPath_ch,genomeName_ch)

    //
    // MODULE: index genome mapped fqs
    //

    samIndex(mappedBam_ch)

    //
    // MODULE: create bedgraph of mapped fqs
    //

    bamCoverage(mappedBam_ch)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
}