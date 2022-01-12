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
  // http://broadinstitute.github.io/picard/explain-flags.html -> -F = remove, 4 = read unmapped, 8 = mate unmapped, 256 = not primary alignment
  if(params.singleEnd) {
    """
    hisat2 -p $task.cpus \
           --very-sensitive \
           --no-spliced-alignment \
           --summary-file ${meta}.hisat2.summary.log \
           -x "${hisat2_indexes}/genome_snp_tran" \
           -U $reads \
            | samtools view -@ ${task.cpus} -bS -F 4 -F 256 - > ${meta}.bam
           
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
            | samtools view -@ ${task.cpus} -bS -F 4 -F 8 -F 256 - > ${meta}.bam
    """
  
  }

}

process SAMTOOLS {
  publishDir "${params.outdir}/${meta}", pattern: '*.sorted.bam.flagstat', mode: 'copy'
  
  tag "$meta"
  
  input:
  tuple val(meta), path(bam)
  
  output:
  path("*.sorted.bam"), emit: bam
  path("*.sorted.bam.flagstat"), emit: flagstat
  
  script:
  """
  samtools sort -@ ${task.cpus} $bam > ${meta}.sorted.bam
  samtools flagstat ${meta}.sorted.bam > ${meta}.sorted.bam.flagstat
  """

}

process PICARD {
  publishDir "${params.outdir}/${sorted_bam.simpleName}", pattern: '*.metrics.txt', mode: 'copy'
  publishDir "${params.outdir}/${sorted_bam.simpleName}", pattern: '*.sorted.md.bam', mode: 'copy'
  publishDir "${params.outdir}/${sorted_bam.simpleName}", pattern: '*.sorted.md.bai', mode: 'copy'
  
  tag "$sorted_bam.simpleName"

  input:
  path(sorted_bam)
  
  output:
  path("*.sorted.md.bam"), emit: bam
  path("*.metrics.txt"), emit: metrics
  tuple path("*.sorted.md.bam"), path("*.sorted.md.bai"), emit: qc
  
  script:
  def avail_mem = 3
  if (!task.memory) {
    log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
  } else {
    avail_mem = task.memory.giga
  }
  """
  picard -Xmx${avail_mem}g \
         MarkDuplicates \
         ASSUME_SORTED=true \
         CREATE_INDEX=true \
         I=$sorted_bam \
         O=${sorted_bam.simpleName}.sorted.md.bam \
         M=${sorted_bam.simpleName}.MarkDuplicates.metrics.txt
  """

}

process FEATURECOUNTS {
  publishDir "${params.outdir}", pattern: '*.txt.gz', mode: 'copy'

  input:
  path(gtf_file_ch)
  path(sorted_bam)
  
  output:
  path("*.txt.gz")
  
  script:
  // ignoreDup = reads that were marked as duplicates will be ignored
  if(params.singleEnd) {
    """
    featureCounts -T $task.cpus \
                  -t exon \
                  -g gene_id \
                  -a $gtf_file_ch \
                  -o featureCounts_output.txt \
		   --ignoreDup \
                  $sorted_bam
                  
    pigz -p $task.cpus featureCounts_output.txt
    """
  } else {
    """
    featureCounts -T $task.cpus \
                  -p \
                  -t exon \
                  -g gene_id \
                  -a $gtf_file_ch \
                  -o featureCounts_output.txt \
                  --ignoreDup \
                  $sorted_bam
                  
    pigz -p $task.cpus featureCounts_output.txt
    """
  
  }
  
}

process DEEPTOOLS {
  publishDir "${params.outdir}/${sorted_bam.simpleName}", pattern: '*.coverage.bw', mode: 'copy'
  
  tag "$sorted_bam.simpleName"
  
  input:
  tuple path(sorted_bam), path(sorted_bam_bai)
  
  output:
  path("*.coverage.bw")
  
  script:
  // generates a coverage track for IGV
  """
  bamCoverage -p $task.cpus \
              -b $sorted_bam \
              --ignoreDuplicates \
              -o ${sorted_bam.simpleName}.coverage.bw
  """
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
  publishDir "${params.outdir}/${sorted_bam.simpleName}/reseq", pattern: "*.{txt,pdf,r,xls}", mode: 'copy'

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
  tuple path(sorted_bam), path(sorted_bam_bai)
    
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
  path(picard_ch)
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

