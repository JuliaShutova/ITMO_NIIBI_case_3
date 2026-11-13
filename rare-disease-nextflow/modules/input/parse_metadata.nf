process PARSE_METADATA {
    tag 'metadata'
    publishDir "${params.outdir}/preprocessing", mode: 'copy'
    input:
    path metadata
    output:
    path 'metadata_parsed.json', emit: parsed
    script:
    """
    python ${projectDir}/scripts/parse_metadata.py ${metadata} metadata_parsed.json
    """
}
