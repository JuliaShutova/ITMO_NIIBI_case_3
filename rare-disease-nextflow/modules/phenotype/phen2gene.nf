process RUN_PHEN2GENE {
    tag 'phen2gene'

    publishDir "${params.outdir}/phenotype/phen2gene", mode: 'copy'

    input:
    path hpo

    output:
    path 'phen2gene_raw/output_file.associated_gene_list', emit: tsv

    script:
    if (!params.phen2gene_script) {
        error "Set params.phen2gene_script to Phen2Gene executable"
    }

    """
    ${params.phen2gene_script} -f ${hpo} -out phen2gene_raw
    """
}

