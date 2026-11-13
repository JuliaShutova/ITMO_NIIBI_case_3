process PREPARE_VARIANTS {
    tag 'prepare-variants'
    publishDir "${params.outdir}/preprocessing", mode: 'copy'
    input:
    path variants
    path metadata
    output:
    path 'variants_prepared.vcf', emit: vcf
    script:
    """
    python ${projectDir}/scripts/prepare_variants.py \
      --input ${variants} --metadata ${metadata} --out variants_prepared.vcf
    """
}
