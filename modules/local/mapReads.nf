#!/usr/bin/env nextflow

params.genomePath='genome_local_placeholder'
params.readsPath='reads_local_placeholder'
params.genomeName='genomeName_local_placeholder'

process mapReads {
	maxForks 3

    publishDir 'output', mode: 'copy'

	input:
		path genomePath
		path readsPath
        val genomeName

	output:
		path "${readsPath.baseName}_mt_${genomeName}.bam"

	script:
	"""
	minimap2 -ax map-ont $genomePath $readsPath --MD |
	samtools view -bS |
	samtools sort -o ${readsPath.baseName}_mt_${genomeName}.bam
	"""
}
