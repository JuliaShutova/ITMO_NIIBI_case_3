process RUN_LIRICAL {
    tag 'lirical'
    publishDir "${params.outdir}/phenotype/lirical", mode: 'copy'
    input:
    path phenopacket
    output:
    path 'lirical.tsv', emit: tsv
    script:
    if (!params.lirical_jar) error "Set params.lirical_jar to LIRICAL.jar"
    if (!params.exomiser_data) error "Set params.exomiser_data to LIRICAL/Exomiser data directory"
    """
    java -jar ${params.lirical_jar} phenopacket \
      --global -e ${params.exomiser_data} -p ${phenopacket} --tsv -x lirical
    test -s lirical.tsv || { echo 'LIRICAL did not create lirical.tsv' >&2; exit 1; }
    """
}
