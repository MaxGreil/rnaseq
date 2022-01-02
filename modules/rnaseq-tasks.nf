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
  publishDir "${params.outdir}/${meta}", pattern: '*.sorted.bam.flagstat', mode: 'copy'
  publishDir "${params.outdir}/${meta}", pattern: '*.sorted.bam.bai', mode: 'copy'

  tag "$meta"
  
  input:
  tuple val(meta), path(mapped_sam)
  
  output:
  val(meta), emit: meta
  path("*.sorted.bam"), emit: bam
  path("*.sorted.bam.flagstat"), emit: flagstat
  tuple path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: all
  
  script:
  """
  samtools view -@ ${task.cpus} -bS -o ${meta}.bam ${mapped_sam}
  samtools sort -@ ${task.cpus} ${meta}.bam > ${meta}.sorted.bam
  samtools flagstat ${meta}.sorted.bam > ${meta}.sorted.bam.flagstat
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

process PRESEQ {
  publishDir "${params.outdir}/${meta}", pattern: '*.txt', mode: 'copy'

  tag "${meta}"
  
  input:
  val(meta)
  tuple path(sorted_bam), path(sorted_bam_bai)
  
  output:
  path("*.txt")
  
  script:
  """
  preseq c_curve -B -o ${meta}.c_curve.txt \
           $sorted_bam
  
  preseq lc_extrap -B -o ${meta}.lc_extrap.txt \
           $sorted_bam
  """

}

process FASTQC {

  tag "${meta}"

  input:
  val(meta)
  path(sorted_bam)
    
  output:
  path('*_fastqc.{zip,html,txt}')
    
  script:
  """
  fastqc -t $task.cpus $sorted_bam
  """

}

process MULTIQC {
  publishDir "${params.outdir}", mode:'copy'

  input:
  path(samtools_ch)
  path(preseq_ch)
  path(fastqc_ch)
  
  output:
  path('multiqc_report.html')

  script:
  """
  multiqc .
  """

}

