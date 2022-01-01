/* 
 * include requires tasks 
 */
include { UNCOMPRESS_GENOTYPE_INDEX; HISAT2; SAMTOOLS; FEATURECOUNTS; } from '../modules/rnaseq-tasks.nf'

/* 
 * define the data analysis workflow 
 */
workflow rnaseqFlow {
    // required inputs
    take:
      reads
    // workflow implementation
    main:
    
      if( params.hisat2_index ) {
        Channel
          .fromPath( params.hisat2_index )
          .ifEmpty { exit 1, "hisat2_index - ${params.hisat2_index} was empty - no input file supplied" }
          .set { hisat2_index_ch }
      }
      
      if( params.gtf_file ) {
        Channel
          .fromPath( params.gtf_file )
          .ifEmpty { exit 1, "gtf_file - ${params.gtf_file} was empty - no input file supplied" }
          .set { gtf_file_ch }
      }
    
      if( params.singleEnd ){
        Channel
          .fromPath( reads )
          .ifEmpty { exit 1, "reads - ${params.reads} was empty - no input files supplied" }
          .map { file -> tuple(file.simpleName, file) }
          .set { reads_ch }
      } else {
        Channel
          .fromFilePairs( reads )
          .ifEmpty { exit 1, "reads - ${params.reads} was empty - no input files supplied" }
          .set { reads_ch }
      }
      
      UNCOMPRESS_GENOTYPE_INDEX(hisat2_index_ch)

      HISAT2(UNCOMPRESS_GENOTYPE_INDEX.out.first(), reads_ch) // value channel, queue channel -> process termination determined by content of queue channel
      
      SAMTOOLS(HISAT2.out)
      
      FEATURECOUNTS(gtf_file_ch, SAMTOOLS.out.bam.collect())
      
      //data <- read.delim("~/Dokumente/rnaseq/output/featureCounts_output.txt", sep = "\t", header=T, skip=1, row.names="Geneid")
      
}
