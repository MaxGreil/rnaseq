/* 
 * include requires tasks 
 */
include { BUILD_HISAT_INDEX; } from '../modules/rnaseq-tasks.nf'

/* 
 * define the data analysis workflow 
 */
workflow rnaseqFlow {
    // required inputs
    take:
      reads
    // workflow implementation
    main:
      if(params.singleEnd){
        Channel
          .fromPath(reads)
          .ifEmpty { exit 1, "${params.reads} was empty - no input files supplied" }
          .set { fastqgz_ch } 
      } else {
        //FilePairs?
      }
      
      fastqgz_ch.view()
      
}
