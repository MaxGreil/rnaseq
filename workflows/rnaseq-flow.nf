/* 
 * include requires tasks 
 */
include { UNCOMPRESS_GENOTYPE_INDEX; HISAT2; SAMTOOLS; UNCOMPRESS_GTF_FILE; STRINGTIE; } from '../modules/rnaseq-tasks.nf'

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
      
      UNCOMPRESS_GTF_FILE(gtf_file_ch)
      
      STRINGTIE(UNCOMPRESS_GTF_FILE.out.first(), SAMTOOLS.out.bam)
      
}
