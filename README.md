# Rare Disease Variant Prioritization Pipeline (`rare-disease-nextflow`)

Автоматизированный биоинформатический пайплайн на базе **Nextflow (DSL2)** для поиска, аннотации и приоритизации генетических вариантов при диагностике редких наследственных заболеваний[cite: 1, 2, 9].  
Решение разработано в рамках **Кейса №3 («Геномный детектив: поиск генетических причин редких наследственных заболеваний»)** от ПИШ ИТМО на ХимБио-Хакатоне «СИНТЕЗ 2026»[cite: 3, 9].

---

## 1. Описание задачи и концепция решения

Секвенирование экзома или генома пациента позволяет выявить тысячи и десятки тысяч генетических вариантов, подавляющее большинство которых являются доброкачественными полиморфизмами или вариантами с неизвестным клиническим значением (VUS)[cite: 9].

Задача системы — сопоставить сырые или предфильтрованные генетические варианты пациента с его клинической картиной (фенотипическими признаками в формате **Human Phenotype Ontology — HPO**), проинтегрировать данные из клинических и популяционных баз и сформировать ранжированный список генов и вариантов-кандидатов с прозрачным объяснением позиции каждого кандидата[cite: 9, 10].

### Ключевые особенности пайплайна:
* **Универсальный вход и строгая валидация:** поддержка входных вариантов в форматах `VCF` и `TSV`, автоматическая проверка целостности входных файлов и нормализация вариантов по референсному геному `GRCh38`[cite: 1, 2, 10].
* **Раздельная обработка классов вариантов:** автоматическая классификация вариантов на малые варианты (`SNV`, `Indel`) и структурные перестройки (`SV`) с маршрутизацией в специализированные ветки аннотации[cite: 1].
* **Комплексная функциональная и клиническая аннотация:** интеграция `SnpEff` / `SnpSift`, локальных баз `ClinVar`, `gnomAD` и `AnnotSV`[cite: 1, 2].
* **Фенотип-ориентированный поиск (Phenotype-driven prioritization):** интеграция инструментов `Phen2Gene` и `LIRICAL` для взвешивания генов по степени соответствия симптомам пациента (HPO)[cite: 1, 2].
* **Отказоустойчивая модульная сборка:** за счет паттерна резервных таблиц (`EMPTY_TABLE`) пайплайн стабильно отрабатывает этап финальной интеграции и скоринга даже в случаях, когда в образце отсутствует определенный класс вариантов (например, нет `SV`) или отключен один из опциональных аннотаторов[cite: 1].
* **Интерпретируемость и работа с неопределенностью:** четкое разделение подтвержденных фактов из баз данных, вычислительных предсказаний и гипотез системы[cite: 10, 11]; поддержка трех базовых клинических сценариев (однозначное моногенное заболевание, конкурирующие гипотезы, недостаточность данных)[cite: 11].

---

## 2. Архитектура пайплайна (`main.nf`)

Воркфлоу построен по стандарту **Nextflow DSL2** и включает следующие модули[cite: 1]:

1. **Входной контроль и парсинг (`modules/input/`):**
   * `VALIDATE_INPUT` — проверка наличия и формата входных файлов (`metadata`, `phenotype`, `variants`)[cite: 1].
   * `PARSE_METADATA` — извлечение клинико-демографических данных и модели наследования[cite: 1, 10].
   * `PARSE_PHENOTYPE` — извлечение HPO-терминов пациента и генерация фенопэкета (`phenopacket`)[cite: 1].
2. **Предобработка и классификация вариантов (`modules/preprocessing/`):**
   * `PREPARE_VARIANTS` — стандартизация входной таблицы/VCF[cite: 1].
   * `NORMALIZE_VCF` — выравнивание и нормализация инделов и мультиаллельных сайтов относительно индексированного референса (`hg38.fa` + `hg38.fa.fai`)[cite: 1, 2].
   * `CLASSIFY_VARIANTS` — разделение потока вариантов на три канала: `snv`, `indel`, `sv` и формирование отчета о классификации[cite: 1].
   * `EXPORT_BASIC_TSV` — экспорт базовых характеристик вариантов по каждому классу[cite: 1].
3. **Функциональная и клиническая аннотация (`modules/annotation/`):**
   * `ANNOTATE_SNPEFF` *(при `enable_snpeff = true`)* — аннотация `SNV` и `Indel` через `SnpEff` / `SnpSift` с подключением базы `ClinVar`[cite: 1, 2].
   * `ANNOTATE_SV` *(при `enable_annotsv = true`)* — аннотация структурных вариантов через `AnnotSV`[cite: 1, 2].
4. **Фенотипический анализ (`modules/phenotype/`):**
   * `RUN_PHEN2GENE` *(при `enable_phen2gene = true`)* — расчет приоритета генов на основе графа HPO-признаков пациента[cite: 1, 2].
   * `RUN_LIRICAL` *(при `enable_lirical = true`)* — оценка правдоподобия кандидатных диагнозов по фенопэкету[cite: 1, 2].
5. **Интеграция, скоринг и отчетность (`modules/integration/` & `modules/reporting/`):**
   * `EMPTY_TABLE` — генерация пустых совместимых структур данных для неактивных веток[cite: 1].
   * `INTEGRATE_AND_RANK` — объединение всех аннотаций, фенотипических весов, метаданных и правил наследования в единую матрицу доказательств и расчет итогового ранга[cite: 1].
   * `GENERATE_PREDICTIONS` — формирование финального ранжированного списка кандидатов с текстовым обоснованием каждого решения[cite: 1, 10].

---

## 3. Структура репозитория

```text
rare-disease-nextflow/
├── main.nf                  # Главный скрипт воркфлоу Nextflow DSL2
├── nextflow.config          # Конфигурация параметров, ресурсов и профилей запуска
├── environment.yml          # Спецификация Conda-окружения
├── modules/                 # Модули пайплайна по этапам обработки
│   ├── input/               # Валидация и парсинг метаданных и HPO
│   ├── preprocessing/       # Подготовка, нормализация VCF и классификация (SNV/Indel/SV)
│   ├── annotation/          # Модули SnpEff и AnnotSV
│   ├── phenotype/           # Модули Phen2Gene и LIRICAL
│   ├── integration/         # Интеграция таблиц и ранжирование кандидатов
│   └── reporting/           # Генерация финальных предсказаний и отчетов
├── data/
│   └── reference/           # Директория для референсного генома (hg38.fa, hg38.fa.fai)
├── inputdata/               # Входные клинические случаи (metadata.json, phenotype.json, variants.tsv/vcf)
└── results/                 # Выходная директория с результатами и отчетами
```

---

## 4. Установка и настройка окружения

### Шаг 1. Создание основного окружения (`rare-disease-nf`)
Создайте и активируйте Conda-окружение из файла `environment.yml`, затем установите `snpsift` и проверьте версию Nextflow:

```bash
conda env create -f environment.yml
conda activate rare-disease-nf
conda install -c bioconda snpsift
nextflow -version
```

### Шаг 2. Создание окружения для Phen2Gene (`phen2gene`)
Инструмент `Phen2Gene` запускается внутри пайплайна в собственном изолированном окружении через `conda run -n phen2gene`[cite: 2]. Создайте его:

```bash
conda env create -f environment.yml -n phen2gene
conda activate phen2gene
# При необходимости доустановите зависимости самого репозитория Phen2Gene
conda activate rare-disease-nf
```

---

## 5. Подготовка референса и баз данных

### 1. Индексация референсного генома (`GRCh38`)
Поместите файл референсного генома `hg38.fa` в папку `data/reference/` (или укажите собственный путь в конфиге или командной строке)[cite: 2]. Пайплайн требует наличия `.fai` индекса рядом с FASTA-файлом[cite: 1]:

```bash
samtools faidx data/reference/hg38.fa
```

### 2. Настройка SnpEff
Для работы модуля функциональной аннотации в `nextflow.config` по умолчанию могут стоять пустые значения (`enable_snpeff = false`, `snpeff_db = ''`, `snpeff_config = ''`). Одного лишь переключения `enable_snpeff = true` недостаточно — должны выполняться **три обязательных условия**:

1. Команда `snpEff` доступна в окружении `rare-disease-nf`.
2. Скачана база данных SnpEff для сборки `GRCh38`.
3. В `nextflow.config` прописаны имя базы (`snpeff_db`) и полный путь к `snpEff.config`[cite: 2].

**Команды для проверки и загрузки базы SnpEff:**
```bash
# Посмотреть список доступных баз для GRCh38
snpEff databases | grep -i GRCh38

# Скачать нужную версию базы (например, GRCh38.105 или GRCh38.99)
snpEff download -v GRCh38.105
```

**Поиск пути к конфигурационному файлу `snpEff.config`:**
```bash
find "$CONDA_PREFIX" -name snpEff.config
```

Скопируйте полученный путь и укажите его в `nextflow.config` в параметре `snpeff_config`, а название скачанной базы — в `snpeff_db`[cite: 2].

---

## 6. Конфигурация (`nextflow.config`)

Пример настроенного блока `params` в файле `nextflow.config`[cite: 2]:

```groovy
params {
    input_dir = 'inputdata'
    metadata  = "${params.input_dir}/metadata.json"
    phenotype = "${params.input_dir}/phenotype.json"
    variants  = "${params.input_dir}/variants.vcf"

    outdir          = 'results'
    assembly        = 'GRCh38'
    reference_fasta = "${projectDir}/data/reference/hg38.fa"
    skip_ref_check  = false

    // Настройки SnpEff и ClinVar
    enable_snpeff = true
    snpeff_db     = 'GRCh38.105' // или 'GRCh38.99'
    snpeff_config = '/path/to/miniconda3/envs/rare-disease-nf/share/snpeff-5.2-3/snpEff.config'
    clinvar_db    = '/path/to/rare-disease-nextflow/clinvar.vcf.gz'

    // Настройки структурных вариантов (AnnotSV)
    enable_annotsv = false
    annotsv_script = ''
    annotsv_db     = ''

    // Настройки Phen2Gene
    enable_phen2gene = true
    phen2gene_script = '/path/to/miniconda3/bin/conda run -n phen2gene python /path/to/rare-disease-nextflow/Phen2Gene/phen2gene.py'
    phen2gene_data   = ''

    // Настройки LIRICAL
    enable_lirical = false
    lirical_jar    = ''
    exomiser_data  = ''

    // Дополнительные локальные базы и ML-модуль
    clinvar_vcf      = ''
    gnomad_vcf       = ''
    gene_disease_tsv = ''
    enable_ml        = false
}
```

### Профили выполнения (`profiles`)[cite: 2]:
* `standard` — локальный запуск (`executor = 'local'`). По умолчанию на задачу выделяется `2 CPU` и `4 GB RAM` (лимит `2h`), а на ресурсоемкие задачи `ANNOTATE_SNPEFF` и `ANNOTATE_SV` — `4 CPU` и `8 GB RAM` (лимит `8h`)[cite: 2].
* `conda` — автоматическое управление зависимостями через `environment.yml`[cite: 2].
* `docker` — запуск внутри Docker-контейнера (`continuumio/miniconda3:24.7.1-0`)[cite: 2].
* `test` — тестовый профиль, отключающий внешние аннотаторы (`enable_snpeff = false`, `enable_annotsv = false`, `enable_phen2gene = false`, `enable_lirical = false`) для быстрой проверки логики пайплайна и парсинга входных данных[cite: 2].

---

## 7. Запуск пайплайна

### Полный боевой запуск на клиническом случае:
Убедитесь, что активно окружение `rare-disease-nf`, и выполните команду:

```bash
conda activate rare-disease-nf

nextflow run main.nf \
  --metadata inputdata/case_01_metadata.json \
  --phenotype inputdata/case_01_phenotype.json \
  --variants inputdata/case_01_variants.tsv \
  --reference_fasta data/reference/hg38.fa \
  --outdir results
```

### Быстрый тестовый прогон (без баз аннотации):
```bash
nextflow run main.nf -profile test \
  --metadata inputdata/case_01_metadata.json \
  --phenotype inputdata/case_01_phenotype.json \
  --variants inputdata/case_01_variants.tsv \
  --reference_fasta data/reference/hg38.fa \
  --outdir results_test
```

### Переопределение параметров из командной строки:
Любой параметр из `nextflow.config` можно передать напрямую при вызове `nextflow run`, например:
```bash
nextflow run main.nf \
  --metadata inputdata/case_02_metadata.json \
  --phenotype inputdata/case_02_phenotype.json \
  --variants inputdata/case_02_variants.vcf \
  --reference_fasta data/reference/hg38.fa \
  --enable_snpeff true \
  --snpeff_db GRCh38.105 \
  --enable_phen2gene true \
  --outdir results/case_02
```

---

## 8. Формат входных и выходных данных

### Входные данные:
1. **`--metadata` (`*.json`)** — метаданные пациента и клинического случая (идентификатор случая, пол, предполагаемый тип наследования, семейный анамнез)[cite: 1, 10].
2. **`--phenotype` (`*.json`)** — перечень клинических симптомов пациента в стандартизированных терминах Human Phenotype Ontology (например, `HP:0001250`, `HP:0001263`)[cite: 1, 10].
3. **`--variants` (`*.tsv` / `*.vcf`)** — таблица или VCF-файл с обнаруженными генетическими вариантами пациента (хромосома, позиция, референсный и альтернативный аллели, генотип, качество)[cite: 1, 2, 10].
4. **`--reference_fasta` (`*.fa` + `*.fa.fai`)** — референсный геном человека сборки `GRCh38` с индексом `samtools`[cite: 1, 2].

### Выходные данные (в директории `--outdir`):
* **Нормализованные и классифицированные варианты:** VCF и базовые TSV-файлы, разделенные по типам (`snv`, `indel`, `sv`), плюс сводный отчет о классификации вариантов[cite: 1].
* **Таблицы аннотаций:** результаты работы `SnpEff` / `ClinVar` (функциональные последствия, затронутые гены и транскрипты, клиническая значимость) и `AnnotSV`[cite: 1, 2, 10].
* **Результаты фенотипического ранжирования:** таблица весов ассоциации «ген — фенотип» от `Phen2Gene` (и `LIRICAL` при включении)[cite: 1].
* **Итоговый рейтинг и предсказания (`INTEGRATE_AND_RANK` / `GENERATE_PREDICTIONS`):**
  * Ранжированный список генов и вариантов-кандидатов с интегральным скором[cite: 1, 10].
  * Детализация вклада каждого фактора (клиническая база `ClinVar`, частота в популяции, тяжесть функционального эффекта `SnpEff`, совпадение с HPO-профилем `Phen2Gene`, соответствие типу наследования)[cite: 1, 2, 10].
  * Интерпретация уровня уверенности и пояснение ограничений (для сложных случаев с конкурирующими гипотезами или неполным фенотипом)[cite: 11].

---

## 9. Ограничения и исследовательский статус

Данное программное обеспечение разработано в рамках хакатона исключительно для **научно-исследовательских и образовательных целей**[cite: 9, 10]. Результаты работы пайплайна (ранги вариантов, оценки патогенности и фенотипического соответствия) не являются медицинским заключением и не могут использоваться для постановки клинического диагноза без верификации врачом-генетиком и лабораторного подтверждения[cite: 10].
