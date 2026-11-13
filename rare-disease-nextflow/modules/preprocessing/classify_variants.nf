process CLASSIFY_VARIANTS {
    tag 'classify'
    publishDir "${params.outdir}/classification", mode: 'copy'
    input:
    path vcf
    output:
    path 'snv.vcf.gz', optional: true, emit: snv
    path 'indel.vcf.gz', optional: true, emit: indel
    path 'sv.vcf.gz', optional: true, emit: sv
    path 'classification_report.tsv', emit: report
    script:
    """
    python ${projectDir}/scripts/classify_variants.py \\
      --vcf ${vcf} --outdir . --report classification_report.tsv
    """
}
