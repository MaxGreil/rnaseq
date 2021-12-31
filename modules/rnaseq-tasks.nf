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
           --downstream-transcriptome-assembly \
           -x "${hisat2_indexes}/genome" \
           -U $reads \
           > ${meta}.sam
    """
  } else {
    """
    hisat2 -p $task.cpus \
           --very-sensitive \
           --no-spliced-alignment \
           --downstream-transcriptome-assembly \
           -x "${hisat2_indexes}/genome" \
           -1 ${reads[0]} \
           -2 ${reads[1]}
           > ${meta}.sam
    """
  }

}

process SAMTOOLS {
  publishDir "${params.outdir}", pattern: "${meta}.sorted.bam", mode: 'copy'
  publishDir "${params.outdir}", pattern: "${meta}.sorted.bam.bai", mode: 'copy'

  tag "$meta"
  
  input:
  tuple val(meta), path(mapped_sam)
  
  output:
  tuple val(meta), path("${meta}.sorted.bam"), emit: bam
  tuple val(meta), path("${meta}.sorted.bam.bai"), emit: bai
  
  script:
  """
  samtools view -@ ${task.cpus} -bS -o ${meta}.bam ${mapped_sam}
  samtools sort -@ ${task.cpus} ${meta}.bam > ${meta}.sorted.bam
  samtools index ${meta}.sorted.bam ${meta}.sorted.bam.bai
  """

}

process UNCOMPRESS_GTF_FILE {

  tag "$gtf_file_ch.simpleName"

  input:
  file(gtf_file_ch)
    
  output:
  path('*')
  
  script:
  """
  gzip -d -f $gtf_file_ch
  """

}

process STRINGTIE {
  publishDir "${params.outdir}", pattern: "${meta}.transcripts.gtf", mode: 'copy'
  
  tag "$meta"
  
  input:
  path(gtf_file_ch)
  tuple val(meta), path(sorted_bam)
  
  output:
  tuple val(meta), path("${meta}.transcripts.gtf")
  
  script:
  """
  stringtie -p $task.cpus \
            -G $gtf_file_ch \
            -o ${meta}.transcripts.gtf \
            -l $meta ${sorted_bam}
  """

}
