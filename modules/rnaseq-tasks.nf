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
  
  tag "$reads.baseName"
  
  input:
  path(hisat2_indexes)
  path(reads)

  output:
  path('*.sam')
  
  script:
  """
  hisat2 -p $task.cpus \
         --very-sensitive \
         --no-spliced-alignment \
         -x "${hisat2_indexes}/genome" \
         -U $reads \
         > ${reads.baseName}.sam
  """

}
