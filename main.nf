#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */
params.outdir = "output"

log.info """\
         reads:          ${params.reads}
         hisat2_indices: ${params.hisat2_indices}
         singleEnd:      ${params.singleEnd}
         outdir:         ${params.outdir}
         tracedir:       ${params.tracedir}
         """
         .stripIndent()

include { rnaseqFlow } from './workflows/rnaseq-flow.nf'

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

    log.info ( workflow.success ? "\nDone! Open the following reports in your browser --> $params.tracedir/execution_report.html\n" : "Oops .. something went wrong" )

}
