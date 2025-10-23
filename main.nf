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
    IMPORT NF-CORE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MINIMAP2_ALIGN } from './modules/nf-core/minimap2/align/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow{

    //
    // CHANNEL: create versions channel
    //

    ch_versions = Channel.empty()

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
        .map { item -> item.fastq_path }
        .flatten()
        .distinct()

    //
    // MODULE: Find number of reads in fastq
    //

    readsCount(
        ch_distinctReads
    )
    ch_fq_count = readsCount.out.count.collect()
    ch_versions = ch_versions.mix(readsCount.out.versions)

    //
    // MODULE: Pool distinct reads into a single csv file
    //

    collectReadCounts(
        ch_fq_count
    )

    //
    // ****************************
    //
    // SECTION: Map all read in fqs to genome
    //
    // ****************************
    //

    //
    // CHANNEL: create input channel
    //

    ch_data = Channel.fromPath(params.inputFile)
                .splitCsv(header: true)
                .map { row ->
                    [[id: row.sample_id, genome_path: row.genome_path], row.fastq_path]
                }

    //
    // MODULE: map fq reads to genome
    //

    MINIMAP2_ALIGN(
        ch_data.map{ meta, fq -> [meta, fq] },
        ch_data.map{ meta, fq -> [meta, meta.genome_path] },
        true,
        false,
        false,
        false
    )
    ch_gen_bam = MINIMAP2_ALIGN.out.bam
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    /*
    mapReads(
        ch_data.map{ meta, fq -> [meta, fq] },
        ch_data.map{ meta, fq -> meta.genome_path }
    )
    ch_gen_bam = mapReads.out.bam
    */

    //
    // MODULE: index genome mapped fqs
    //

    samIndex(
        ch_gen_bam.map{ meta, bam -> [meta, bam] }
    )
    ch_gen_bai = samIndex.out.bai

    //
    // CHANNEL: Combine BAM and BAI
    //
    ch_gen_bam_bai = ch_gen_bam
        .join(ch_gen_bai, by: [0])
        .map {
            meta, bam, bai ->
                if (bai) {
                    [ meta, bam, bai ]
                }
        }

    //
    // CHANNEL: Filter empty bams
    //
    ch_gen_bam_bai = ch_gen_bam_bai.filter { row -> 
            file(row[1]).size() >= params.min_bam_size 
            }

    //
    // MODULE: create bedgraph of mapped fqs
    //

    bamCoverage(
        ch_gen_bam_bai.map{ meta, bam, bai -> [meta, bam, bai] }
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
}