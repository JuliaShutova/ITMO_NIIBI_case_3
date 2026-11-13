process EXPORT_BASIC_TSV {
    tag "basic-${category}"
    publishDir "${params.outdir}/preprocessing/variant_tables", mode: 'copy'
    input:
    tuple val(category), path(vcf)
    output:
    path "${category}_basic.tsv", emit: tsv
    script:
    """
    python ${projectDir}/scripts/vcf_to_basic_tsv.py --vcf ${vcf} --out ${category}_basic.tsv
    """
}
