process VALIDATE_INPUT {
    tag 'input-qc'
    publishDir "${params.outdir}/qc", mode: 'copy'
    input:
    path metadata
    path phenotype
    path variants
    output:
    path 'input_qc_report.json'
    script:
    """
    python ${projectDir}/scripts/validate_input.py \
      --metadata ${metadata} --phenotype ${phenotype} --variants ${variants} \
      --out input_qc_report.json
    """
}
