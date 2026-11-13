nextflow.enable.dsl=2

include { VALIDATE_INPUT } from './modules/input/validate_input'
include { PARSE_METADATA } from './modules/input/parse_metadata'
include { PARSE_PHENOTYPE } from './modules/input/parse_phenotype'
include { PREPARE_VARIANTS } from './modules/preprocessing/prepare_variants'
include { NORMALIZE_VCF } from './modules/preprocessing/normalize_vcf'
include { CLASSIFY_VARIANTS } from './modules/preprocessing/classify_variants'
include { EXPORT_BASIC_TSV } from './modules/preprocessing/export_basic_tsv'
include { ANNOTATE_SNPEFF } from './modules/annotation/snpeff'
include { ANNOTATE_SV } from './modules/annotation/annotsv'
include { RUN_PHEN2GENE } from './modules/phenotype/phen2gene'
include { RUN_LIRICAL } from './modules/phenotype/lirical'
include { INTEGRATE_AND_RANK } from './modules/integration/integrate_and_rank'
include { EMPTY_TABLE } from './modules/integration/empty_table'
include { GENERATE_PREDICTIONS } from './modules/reporting/predictions'

workflow {
    metadata = file(params.metadata, checkIfExists: true)
    phenotype = file(params.phenotype, checkIfExists: true)
    variants = file(params.variants, checkIfExists: true)

    VALIDATE_INPUT(metadata, phenotype, variants)
    PARSE_METADATA(metadata)
    PARSE_PHENOTYPE(phenotype)
    PREPARE_VARIANTS(variants, metadata)

    reference_fasta = file(params.reference_fasta, checkIfExists: true)
    reference_fai = file("${params.reference_fasta}.fai", checkIfExists: true)

    NORMALIZE_VCF(
    PREPARE_VARIANTS.out.vcf,
    reference_fasta,
    reference_fai
    )
    
    CLASSIFY_VARIANTS(NORMALIZE_VCF.out.vcf)

    // Optional outputs are backed by an explicit empty TSV so integration can run
    // even when a variant class or optional annotation tool is absent.
    EMPTY_TABLE()

    snv_indel = CLASSIFY_VARIANTS.out.snv.map { tuple('snv', it) }
        .mix(CLASSIFY_VARIANTS.out.indel.map { tuple('indel', it) })
    sv = CLASSIFY_VARIANTS.out.sv.map { tuple('sv', it) }
    classified_all = snv_indel.mix(sv)
    EXPORT_BASIC_TSV(classified_all)

    if (params.enable_snpeff) {
        ANNOTATE_SNPEFF(snv_indel)
        small_tsv = ANNOTATE_SNPEFF.out.tsv.mix(EMPTY_TABLE.out.small)
    } else {
        small_tsv = EMPTY_TABLE.out.small
    }

    if (params.enable_annotsv) {
        ANNOTATE_SV(sv)
        sv_tsv = ANNOTATE_SV.out.tsv.mix(EMPTY_TABLE.out.sv)
    } else {
        sv_tsv = EMPTY_TABLE.out.sv
    }

    if (params.enable_phen2gene) {
    RUN_PHEN2GENE(PARSE_PHENOTYPE.out.present_hpo)
    p2g_tsv = RUN_PHEN2GENE.out.tsv
    } else {
    p2g_tsv = EMPTY_TABLE.out.phen2gene
    }

    if (params.enable_lirical) {
        RUN_LIRICAL(PARSE_PHENOTYPE.out.phenopacket)
        lirical_tsv = RUN_LIRICAL.out.tsv.mix(EMPTY_TABLE.out.lirical)
    } else {
        lirical_tsv = EMPTY_TABLE.out.lirical
    }

    INTEGRATE_AND_RANK(
        small_tsv.collect(),
        sv_tsv.collect(),
        p2g_tsv.collect(),
        lirical_tsv.collect(),
        EXPORT_BASIC_TSV.out.tsv.collect(),
        CLASSIFY_VARIANTS.out.report,
        PARSE_METADATA.out.parsed,
        PARSE_PHENOTYPE.out.parsed
    )

    GENERATE_PREDICTIONS(INTEGRATE_AND_RANK.out.ranked, PARSE_METADATA.out.parsed)
}
