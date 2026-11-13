
process ANNOTATE_SNPEFF {
    tag { category }

    publishDir { "${params.outdir}/annotation/${category}" }, mode: 'copy'

    input:
    tuple val(category), path(vcf)

    output:
    path "${category}_annotated.vcf", emit: vcf
    path "${category}_annotated.tsv", emit: tsv

    script:
    if (!params.snpeff_db) {
        error "params.snpeff_db must be set when enable_snpeff=true"
    }

    """
    snpEff -Xmx8g \
        -c ${params.snpeff_config ?: '$HOME/.config/snpEff/snpEff.config'} \
        -v ${params.snpeff_db} \
        ${vcf} > ${category}_annotated.vcf

    SnpSift annotate \
    -info CLNSIG,CLNDN,CLNREVSTAT,CLNVC,ALLELEID \
    -name CLINVAR_ \
    "${params.clinvar_db}" \
    "${category}_annotated.vcf" \
    > "${category}_clinvar.vcf"


    SnpSift extractFields \
    "${category}_clinvar.vcf" \
    CHROM POS ID REF ALT QUAL FILTER \
    "ANN[*].EFFECT" \
    "ANN[*].IMPACT" \
    "ANN[*].GENE" \
    "ANN[*].GENEID" \
    "ANN[*].FEATURE" \
    "ANN[*].BIOTYPE" \
    "ANN[*].HGVS_C" \
    "ANN[*].HGVS_P" \
    "CLINVAR_CLNSIG" \
    "CLINVAR_CLNDN" \
    "CLINVAR_CLNREVSTAT" \
    "CLINVAR_CLNVC" \
    "CLINVAR_ALLELEID" \
    > "${category}_annotated.tsv"

    """
}
