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

process HISAT2_TO_BAM {
  publishDir "${params.outdir}/${meta}", pattern: '*.log', mode: 'copy'
  
  tag "$meta"
  
  input:
  path(hisat2_indexes)
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path("*.bam"), emit: bam
  path("*.log"), emit: log
  
  script:
  if(params.singleEnd) {
    """
    hisat2 -p $task.cpus \
           --very-sensitive \
           --no-spliced-alignment \
           --summary-file ${meta}.hisat2.summary.log \
           -x "${hisat2_indexes}/genome_snp_tran" \
           -U $reads \
            | samtools view -@ ${task.cpus} -bS - > ${meta}.bam
           
    """
  } else {
    """
    hisat2 -p $task.cpus \
           --very-sensitive \
           --no-spliced-alignment \
           --summary-file ${meta}.hisat2.summary.log \
           -x "${hisat2_indexes}/genome_snp_tran" \
           -1 ${reads[0]} \
           -2 ${reads[1]} \
            | samtools view -@ ${task.cpus} -bS - > ${meta}.bam
    """
  }

}

process SAMTOOLS {
  publishDir "${params.outdir}/${meta}", pattern: '*.sorted.bam', mode: 'copy'
  publishDir "${params.outdir}/${meta}", pattern: '*.sorted.bam.flagstat', mode: 'copy'
  publishDir "${params.outdir}/${meta}", pattern: '*.sorted.bam.bai', mode: 'copy'

  tag "$meta"
  
  input:
  tuple val(meta), path(bam)
  
  output:
  path("*.sorted.bam"), emit: bam
  path("*.sorted.bam.flagstat"), emit: flagstat
  tuple path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: qc
  
  script:
  """
  samtools sort -@ ${task.cpus} $bam > ${meta}.sorted.bam
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
  publishDir "${params.outdir}/${sorted_bam.simpleName}", pattern: '*.txt', mode: 'copy'

  tag "$sorted_bam.simpleName"
  
  input:
  tuple path(sorted_bam), path(sorted_bam_bai)
  
  output:
  path("*.txt")
  
  script:
  """
  preseq c_curve -B -o ${sorted_bam.simpleName}.c_curve.txt \
           $sorted_bam
  
  preseq lc_extrap -B -o ${sorted_bam.simpleName}.lc_extrap.txt \
           $sorted_bam
  """

}

process UNCOMPRESS_BED {

  tag "$bed_file_ch.simpleName"

  input:
  path(bed_file_ch)
  
  output:
  path('*')
  
  script:
  """
  gunzip -d -f $bed_file_ch
  """

}

process RSEQC {

  tag "$sorted_bam.simpleName"
  
  input:
  path(bed_file_ch)
  tuple path(sorted_bam), path(sorted_bam_bai)
  
  output:
  path("*.{txt,pdf,r,xls}")
  
  script:
  """

  read_distribution.py -i $sorted_bam \
                       -r $bed_file_ch \
                        > ${sorted_bam.simpleName}.read_dist.txt

  read_duplication.py -i $sorted_bam \
                      -o ${sorted_bam.simpleName}.read_duplication

  infer_experiment.py -i $sorted_bam \
                      -r $bed_file_ch \
                       > ${sorted_bam.simpleName}.infer_experiment.txt
  """

}

process FASTQC {
  
  tag "$sorted_bam.simpleName"
  
  input:
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
  path(hisat2_ch)
  path(fastqc_ch)
  path(samtools_ch)
  path(preseq_ch)
  path(rseqc_ch)

  output:
  path('multiqc_report.html')

  script:
  """
  multiqc .
  """

}

