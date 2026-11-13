process EMPTY_TABLE {
    output:
    path 'empty_small.tsv', emit: small
    path 'empty_sv.tsv', emit: sv
    path 'empty_phen2gene.tsv', emit: phen2gene
    path 'empty_lirical.tsv', emit: lirical
    script:
    """
    printf 'variant_id\\tgene\\tconsequence\\timpact\\n' > empty_small.tsv
    cp empty_small.tsv empty_sv.tsv
    cp empty_small.tsv empty_phen2gene.tsv
    cp empty_small.tsv empty_lirical.tsv
    """
}
