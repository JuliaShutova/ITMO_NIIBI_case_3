process ANNOTATE_SV {
    tag 'sv'
    publishDir "${params.outdir}/annotation/sv", mode: 'copy'
    input:
    tuple val(category), path(vcf)
    output:
    path 'sv_annotated.tsv', emit: tsv
    script:
    if (!params.annotsv_script) error "Set params.annotsv_script to the AnnotSV executable/script"
    """
    ${params.annotsv_script} -SVinputFile ${vcf} -outputFile sv_annotated.tsv
    """
}
