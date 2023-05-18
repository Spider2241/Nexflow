nextflow.enable.dsl = 2

params.accession = "null"
params.outdir = "SRA_data"


process sra_prefetch {
  storeDir "${params.outdir}"
  container "https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0"

  input:
    val accession

  output:
    path "sra/${accession}.sra", emit: srafile 

  script:
    """
    prefetch $accession 
    """
}

process fastqconvert {
  storeDir "${params.outdir}/sra/${accession}_data/raw_fastq"

  container "https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0"

  input: 
    val accession
    path srafile

  output:
     path '*.fastq', emit: fastqfiles

  script:
  """
  fastq-dump --split-files ${srafile}
  """
}


process fastp {
  publishDir "${params.outdir}/sra/${accession}_data/", mode: 'copy', overwrite: true
  
  container "https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0"

  input:
    path fastqfiles
    val accession

  output:
    path "fastp_fastq/*.fastq", emit: fastqfiles
    path "fastp_report", emit: fastpreport

  script:
  
      if(fastqfiles instanceof List) {
      """
      mkdir fastp_fastq
      mkdir fastp_report
      fastp -i ${fastqfiles[0]} -I ${fastqfiles[1]} -o fastp_fastq/${fastqfiles[0].getSimpleName()}_fastp.fastq -O fastp_fastq/${fastqfiles[1].getSimpleName()}_fastp.fastq -h fastp_report/fastp.html -j fastp_report/fastp.json
      """
    } else {
      """
      mkdir fastp_fastq
      mkdir fastp_report
      fastp -i ${fastqfiles} -o fastp_fastq/${fastqfiles.getSimpleName()}_fastp.fastq -h fastp_report/fastp.html -j fastp_report/fastp.json
      """
    }
}

process fastqc {
  publishDir "${params.outdir}/sra/${params.accession}_data/", mode: 'copy', overwrite: true

  container "https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0"

  input:
    path fastqfiles

  output:
    path "fastqc_report/*.html", emit: results
    path "fastqc_report/*.zip", emit: fastqczip

  script:
    """
    mkdir fastqc_report
    fastqc ${fastqfiles} --outdir fastqc_report
    """
}



workflow {
  srafile = sra_prefetch(params.accession)
  converted = fastqconvert(params.accession, srafile)
  trimmed = fastp(converted, params.accession).fastqfiles
  trimmed_flat = trimmed.flatten()
  report = fastqc(trimmed_flat.collect()).results
}


