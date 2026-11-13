process GENERATE_PREDICTIONS {
    tag 'predictions'
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path ranked
    path metadata
    output:
    path 'predictions.tsv'
    script:
    """
    python ${projectDir}/scripts/predictions.py \
      --ranked ${ranked} --metadata ${metadata} --out predictions.tsv
    """
}
