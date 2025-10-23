#!/usr/bin/env nextflow

process mapReads {

	input:
		tuple val(meta), path(readsPath)
		path genomePath

	output:
		tuple val(meta), path("${readsPath.baseName}_mt_${genomePath.baseName}.bam"), emit: bam

	script:
	"""
	minimap2 -ax map-ont $genomePath $readsPath --MD |
	samtools view -bS |
	samtools sort -o ${readsPath.baseName}_mt_${genomePath.baseName}.bam
	"""
}
