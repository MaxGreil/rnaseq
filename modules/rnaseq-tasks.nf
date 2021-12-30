process UNCOMPRESS_GENOTYPE_INDEX {

  tag "$hisat2_index_ch.baseName"

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
           -x "${hisat2_indexes}/genome" \
           -U $reads \
           > ${meta}.sam
    """
  } else {
    """
    hisat2 -p $task.cpus \
           --very-sensitive \
           --no-spliced-alignment \
           -x "${hisat2_indexes}/genome" \
           -1 ${reads[0]} \
           -2 ${reads[1]}
           > ${meta}.sam
    """
  }

}

process SAMTOOLS {
  publishDir "${params.outdir}", pattern: "${meta}.sorted.bam", mode: 'copy'

  tag "$meta"
  
  input:
  tuple val(meta), path(mapped_sam)
  
  output:
  tuple val(meta), path("${meta}.sorted.bam")
  
  script:
  """
  samtools view -@ ${task.cpus} -bS -o ${meta}.bam ${mapped_sam}
  samtools sort -@ ${task.cpus} ${meta}.bam > ${meta}.sorted.bam
  """

}
