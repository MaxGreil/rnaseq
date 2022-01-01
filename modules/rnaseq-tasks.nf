process UNCOMPRESS_GENOTYPE_INDEX {

  tag "$hisat2_index_ch.simpleName"

  input:
  path(hisat2_index_ch)
  
  output:
  path('*')
  
  script:
  """
  tar xvzf $hisat2_index_ch
  """

}

process HISAT2 {
  
  tag "$meta"
  
  input:
  path(hisat2_indexes)
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path('*.sam')
  
  script:
  if(params.singleEnd) {
    """
    hisat2 -p $task.cpus \
           --very-sensitive \
           --no-spliced-alignment \
           -x "${hisat2_indexes}/genome_snp_tran" \
           -U $reads \
           > ${meta}.sam
    """
  } else {
    """
    hisat2 -p $task.cpus \
           --very-sensitive \
           --no-spliced-alignment \
           -x "${hisat2_indexes}/genome_snp_tran" \
           -1 ${reads[0]} \
           -2 ${reads[1]}
           > ${meta}.sam
    """
  }

}

process SAMTOOLS {
  publishDir "${params.outdir}/${meta}", pattern: '*.sorted.bam', mode: 'copy'
  publishDir "${params.outdir}/${meta}", pattern: '*.sorted.bam.bai', mode: 'copy'

  tag "$meta"
  
  input:
  tuple val(meta), path(mapped_sam)
  
  output:
  path("*.sorted.bam"), emit: bam
  path("*.sorted.bam.bai"), emit: bai
  
  script:
  """
  samtools view -@ ${task.cpus} -bS -o ${meta}.bam ${mapped_sam}
  samtools sort -@ ${task.cpus} ${meta}.bam > ${meta}.sorted.bam
  samtools index ${meta}.sorted.bam ${meta}.sorted.bam.bai
  """

}

process FEATURECOUNTS {
  publishDir "${params.outdir}", pattern: '*.txt', mode: 'copy'

  input:
  path(gtf_file_ch)
  path(sorted_bam)
  
  output:
  path("*.txt")
  
  script:
  if(params.singleEnd) {
    """
    featureCounts -T $task.cpus \
                  -t exon \
                  -g gene_id \
                  -a $gtf_file_ch \
                  -o featureCounts_output.txt \
                  $sorted_bam
    """
  } else {
    """
    featureCounts -T $task.cpus \
                  -p \
                  -t exon \
                  -g gene_id \
                  -a $gtf_file_ch \
                  -o featureCounts_output.txt \
                  $sorted_bam
    """
  
  }
}
