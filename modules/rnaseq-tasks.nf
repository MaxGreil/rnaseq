process HISAT {
  
  tag "$reads.baseName"
  
  input:
  path(reads)

  output:
  stdout
  
  script:
  """
  echo ${task.memory} ${task.cpus}
  """

}
