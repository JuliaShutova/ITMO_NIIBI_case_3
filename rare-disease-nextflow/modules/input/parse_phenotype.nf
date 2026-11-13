process PARSE_PHENOTYPE {
    tag 'phenotype'
    publishDir "${params.outdir}/preprocessing", mode: 'copy'
    input:
    path phenotype
    output:
    path 'phenotype_parsed.json', emit: parsed
    path 'present_hpo.txt', emit: present_hpo
    path 'excluded_hpo.txt', emit: excluded_hpo
    path 'phenopacket.json', emit: phenopacket
    script:
    """
    python ${projectDir}/scripts/parse_phenotype.py ${phenotype} \
      --parsed phenotype_parsed.json --present present_hpo.txt \
      --excluded excluded_hpo.txt --phenopacket phenopacket.json
    """
}
