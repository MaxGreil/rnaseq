/* 
 * include requires tasks 
 */
include { HISAT; } from '../modules/rnaseq-tasks.nf'

/* 
 * define the data analysis workflow 
 */
workflow rnaseqFlow {
    // required inputs
    take:
      reads
    // workflow implementation
    main:
      if( params.singleEnd ){
        Channel
          .fromPath( reads )
          .ifEmpty { exit 1, "${params.reads} was empty - no input files supplied" }
          .set { reads_ch } 
      } else {
        Channel
          .fromFilePairs( reads )
          .ifEmpty { exit 1, "${params.reads} was empty - no input files supplied" }
          .set { reads_ch } 
      }
      
      HISAT(reads_ch)
      
      HISAT.out.view()
      
}
