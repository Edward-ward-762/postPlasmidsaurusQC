#!/usr/bin/env nextflow

process mapReads {

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
