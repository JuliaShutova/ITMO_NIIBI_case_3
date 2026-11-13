process INTEGRATE_AND_RANK {
    tag 'integrate-rank'
    publishDir "${params.outdir}/prioritization", mode: 'copy'
    publishDir "${params.outdir}/integration", mode: 'copy', pattern: 'integrated_*.tsv'
    input:
    path small_tables
    path sv_tables
    path p2g_tables
    path lirical_tables
    path base_tables
    path classification_report
    path metadata
    path phenotype
    output:
    path 'ranked_variants.tsv', emit: ranked
    path 'integrated_variants.tsv'
    path 'integration_report.json'
    script:
    """
    python ${projectDir}/scripts/integrate_and_rank.py \
      --small "${small_tables}" --sv "${sv_tables}" --phen2gene "${p2g_tables}" \
      --lirical "${lirical_tables}" --base "${base_tables}" --classification ${classification_report} \
      --metadata ${metadata} --phenotype ${phenotype} \
      --out ranked_variants.tsv --integrated integrated_variants.tsv \
      --report integration_report.json
    """
}
