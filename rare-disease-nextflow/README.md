# Rare-disease Nextflow pipeline

Исследовательский pipeline для первичной обработки вариантов и evidence-based приоритизации. Не является клинической диагностической системой.

## 1. Требования

- Linux x86_64 (рекомендуется Ubuntu 22.04/24.04)
- Java 17+
- Nextflow 24.10+ (DSL2)
- Conda/Mamba
- Для VCF-обработки: bcftools, bgzip, tabix
- Для реальной аннотации: локальная база SnpEff; для SV — AnnotSV и его ресурсы
- Для фенотипических инструментов: установленные Phen2Gene и LIRICAL плюс их базы

## 2. Входы

Положите в `inputdata/`:
- `metadata.json`
- `phenotype.json`
- `variants.vcf` (стандартный VCF **или** TSV с колонками `chrom,pos,ref,alt,genotype`)

TSV будет преобразован в VCF. Пример TSV не является VCF.

## 3. Установка

```bash
conda env create -f environment.yml
conda activate rare-disease-nf
nextflow -version
bcftools --version
```

Установите Nextflow по официальной инструкции Nextflow. Подготовьте базы отдельно и укажите пути в `nextflow.config` или CLI-параметрах.

## 4. Тест

```bash
python -m unittest discover -s tests -v
nextflow run main.nf -profile test
```

Профиль `test` проверяет ядро pipeline без внешних баз аннотации. Он не имитирует ClinVar, gnomAD, Phen2Gene или LIRICAL.

## 5. Полный запуск

Сначала настройте в `nextflow.config`:
- `enable_snpeff = true`, `snpeff_db = 'имя_установленной_базы'`
- `enable_annotsv = true`, `annotsv_script = '/path/to/AnnotSV'`
- `enable_phen2gene = true`, `phen2gene_script = '/path/to/phen2gene.py'`
- `enable_lirical = true`, `lirical_jar = '/path/to/LIRICAL.jar'`, `exomiser_data = '/path/to/data'`

Затем:
```bash
nextflow run main.nf -profile standard
```

Параметры можно переопределять CLI:
```bash
nextflow run main.nf -profile standard \
  --input_dir inputdata --outdir results \
  --enable_snpeff true --snpeff_db GRCh38.105
```

## 6. Результаты

- `results/qc/input_qc_report.json`
- `results/preprocessing/metadata_parsed.json`
- `results/preprocessing/phenotype_parsed.json`
- `results/preprocessing/variants_normalized.vcf.gz`
- `results/classification/`
- `results/preprocessing/variant_tables/` — базовые строки для всех найденных вариантов, даже если внешняя аннотация выключена
- `results/annotation/` — только если включены соответствующие аннотаторы
- `results/phenotype/` — только если включены Phen2Gene/LIRICAL
- `results/prioritization/ranked_variants.tsv`
- `results/predictions.tsv`

## 7. Важные ограничения текущей версии

1. SnpEff process выполняет функциональную аннотацию и выгружает стандартные ANN-поля. ClinVar и gnomAD не подключаются автоматически в этой базовой версии: локальные VCF этих ресурсов требуют отдельного `bcftools annotate` с корректно проверенными INFO tags/headers и согласованными версиями. Не считать эти поля аннотированными, если они не появились в VCF.
2. OMIM требует соблюдения условий доступа/лицензии. Для gene-disease интеграции в этой версии предусмотрен дальнейший extension point, но автоматическая загрузка OMIM не реализована.
3. Phen2Gene output format зависит от релиза; нормализатор сохраняет табличный output, но его столбцы нужно проверить на установленной версии.
4. Phenopacket, создаваемый парсером, минимальный. Проверьте его соответствие версии Phenopacket Schema, которую принимает ваша сборка LIRICAL. Не передавайте LIRICAL этот JSON без валидации.
5. LIRICAL TSV может называться иначе в разных релизах. Процесс проверяет `lirical.tsv`; при несовпадении имени адаптируйте `-x`/имя выходного файла по `java -jar LIRICAL.jar --help`.
6. В этой версии Phen2Gene и LIRICAL запускаются независимо; их результаты не смешиваются как независимые вероятности. Скрипт baseline использует Phen2Gene gene score только при распознавании соответствующей колонки.
7. Rule-based score — эвристика, не вероятность патогенности. На шести development cases нельзя обучить надёжную универсальную ML-модель.
8. Текущий код не реализует клинически валидированную интеграцию ClinVar/gnomAD/OMIM и не должен использоваться для клинических решений.

## 8. Что следует сделать перед production/benchmark

- Закрепить версии всех инструментов и баз.
- Встроить ClinVar/gnomAD annotation с проверкой build, contig naming и INFO headers.
- Проверить LIRICAL/Phen2Gene CLI на конкретных релизах.
- Подключить gene-disease association source с provenance.
- Добавить полноценную схему Phenopacket и валидацию.
- Провести benchmark на development cases без утечки в hidden cases.
- Проверить, что variant key нормализован одинаково во всех таблицах.
