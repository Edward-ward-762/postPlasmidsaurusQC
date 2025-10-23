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

include { readsCount             } from './modules/local/readsCount.nf'
include { collectReadCounts      } from './modules/local/collectReadCounts.nf'
include { DUMP_SOFTWARE_VERSIONS } from './modules/local/dump_software_versions.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MINIMAP2_ALIGN        } from './modules/nf-core/minimap2/align/main.nf'
include { SAMTOOLS_INDEX        } from './modules/nf-core/samtools/index/main.nf'
include { DEEPTOOLS_BAMCOVERAGE } from './modules/nf-core/deeptools/bamcoverage/main.nf'

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

    //
    // MODULE: index genome mapped fqs
    //

    SAMTOOLS_INDEX(
        ch_gen_bam.map{ meta, bam -> [meta, bam] }
    )
    ch_gen_bai = SAMTOOLS_INDEX.out.bai
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

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

    DEEPTOOLS_BAMCOVERAGE(
        ch_gen_bam_bai.map{ meta, bam, bai -> [meta, bam, bai] },
        [],
        [],
        [[],[]]
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions)

    //
    // ****************************
    //
    // SECTION: software version reporting
    //
    // ****************************
    //

    //
    // MODULE: Collect software versions
    //

    /*
    DUMP_SOFTWARE_VERSIONS (
        ch_versions.unique().collectFile()
    )
    */
    
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
}