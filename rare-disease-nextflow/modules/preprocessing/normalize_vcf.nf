process NORMALIZE_VCF {

    tag 'normalize'

    publishDir "${params.outdir}/preprocessing", mode: 'copy'

    input:
    path vcf
    path reference_fasta
    path reference_fai

    output:
    path 'variants_normalized.vcf.gz', emit: vcf
    path 'variants_normalized.vcf.gz.tbi', emit: index

    script:
    """
    # Добавляем определения contig из индекса референса
    bcftools reheader \
        -f ${reference_fai} \
        -o variants_reheadered.vcf \
        ${vcf}

    # Нормализуем варианты относительно GRCh38
    bcftools norm \
        -f ${reference_fasta} \
        -m -any \
        variants_reheadered.vcf \
        -Oz \
        -o variants_normalized.vcf.gz

    # Создаём индекс
    bcftools index -t variants_normalized.vcf.gz
    """
}
