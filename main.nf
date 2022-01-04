#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */
params.reads = "$baseDir/data/*.fastq.gz" //default: "$baseDir/data/*.fastq.gz" + params.singleEnd = true, alternative: "$baseDir/data/*_{1,2}*.fastq.gz" + params.singleEnd = false
params.outdir = "output"

log.info """\
         reads:        ${params.reads}
         hisat2_index: ${params.hisat2_index}
         gtf_file:     ${params.gtf_file}
         bed_file:     ${params.bed_file}
         singleEnd:    ${params.singleEnd}
         outdir:       ${params.outdir}
         tracedir:     ${params.tracedir}
         """
         .stripIndent()

include { rnaseqFlow } from './workflows/rnaseq-flow.nf'

if( !params.reads ) {

  exit 1, "\nPlease give in reads via --reads <location> \n"

}

/*
 * main script flow
 */
workflow {
	
    rnaseqFlow( params.reads )

}

/*
 * completion handler
 */
workflow.onComplete {

    log.info ( workflow.success ? "\nDone! Open the following reports in your browser --> $params.outdir/multiqc_report.html and $params.tracedir/execution_report.html\n" : "Oops .. something went wrong" )

}
