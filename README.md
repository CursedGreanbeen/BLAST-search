/scripts/BLAST-search/
/computer_data/scripts/

# BLAST-search

Локальный поиск BLASTP для белковых комплексов в формате FASTA (допускаются многоцепочечные файлы). Поиск идёт по локально загруженной базе PDB (pdbaa). ANARCI определяет CDR-H3 в тяжёлых цепях антител. Попадания фильтруются по проценту идентичности внутри CDR-H3, а не по идентичности всей цепи. Результаты агрегируются по уникальным последовательностям CDR-H3, то есть для каждого H3 видно, какие PDB ID ему соответствуют.

## Структура / рабочий процесс

Небольшой конвейер с одной точкой входа:

```
<fasta_dir>/*.fasta
   └─► run_cdr_h3_search.py          (CLI-точка входа / оркестратор)
         ├─► blast_pdb_search_local.py  (библиотека: разбор FASTA и дедупликация,
         │                               поиск CDR-H3 через ANARCI, воркеры blastp,
         │                               разбор XML, фильтрация по идентичности —
         │                               без __main__, напрямую не запускается)
         └─► output_dir/  *.xml, *_filtered_hits.csv, combined_filtered_hits.csv,
                          unique_cdr_h3_summary.csv, errors.txt
```

Скрипты выше не импортируют `normalize_cdr_h3_length.py`. Это отдельный вспомогательный модуль. Он проверяет, одинаковы ли два H3. Его используют `../analyse_cdr_h3/test_h3_threshold.py` и соседний пакет `../blast-search-levenshtein/` (там лежит точная копия). Левенштейн-пакет использует фильтр не по статичному проценту идентичности, а по расстоянию Левенштейна с учётом длины (через rapidfuzz). Далее `../analyse_cdr_h3/analyze_hits.py` читает файл `combined_filtered_hits.csv` этого пакета.

## Требования

- В `PATH` должны быть `blastp` и `blastdbcmd` (NCBI BLAST+).
- Нужна локальная база pdbaa: `update_blastdb.pl --decompress pdbaa`. Если чего-то из этого нет, скрипт сразу завершается с ошибкой.
- Пакеты Python: Biopython и anarci (ANARCI).
- В проектном интерпретаторе `/opt/conda/bin/python3.10` есть Biopython и rapidfuzz, но нет anarci. Для запуска нужно окружение с anarci, например `conda install -c bioconda anarci`.

## Запуск

### run_cdr_h3_search.py — единственный запускаемый скрипт

Запускайте из каталога пакета. Модуль `blast_pdb_search_local` импортируется как локальный, а пути вывода по умолчанию относительные.

```
cd /mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/scripts/BLAST-search
/opt/conda/bin/python3.10 run_cdr_h3_search.py <fasta_dir> \
    --db /path/to/db/pdbaa \
    --workers 4 --threads-per-search 2 \
    --min-h3-identity 80 --max-h3-identity 100 \
    --evalue 10 \
    --output-dir results/blast_pdb_search_local
```

### Параметры (как заданы в argparse)

| Параметр | По умолчанию | Назначение |
|---|---|---|
| `fasta_dir` (позиционный) | — | Каталог с файлами `.fasta`/`.fa`/`.faa`; многозаписные комплексы допускаются |
| `--db` | обязательный | Путь к локальной базе BLAST (например, `.../pdbaa`) |
| `--workers` | 4 | Число одновременно работающих процессов blastp (workers × threads ≤ числа ядер) |
| `--threads-per-search` | 2 | Значение `-num_threads` для каждого blastp |
| `--min-identity` / `--max-identity` | 0 / 100 | Границы идентичности всей цепи в %. Значения выводятся в отчёт, но для запросов с H3 больше не участвуют в фильтрации |
| `--min-h3-identity` / `--max-h3-identity` | 80 / 100 | Границы идентичности CDR-H3 в % — именно по ним идёт фильтрация |
| `--require-h3` | выкл. | Отбрасывать последовательности, в которых ANARCI не нашёл CDR-H3. Параметр разбирается, но в очередь BLAST и так попадают только последовательности с H3, поэтому на поиск он дальше не влияет |
| `--evalue` | 10.0 | E-value для blastp |
| `--output-dir` | `results/blast_pdb_search_local` | Каталог для результатов |

### Порядок работы

1. Разобрать все цепи и убрать дубликаты по SHA1.
2. Параллельно прогнать ANARCI, чтобы найти CDR-H3 в тяжёлых цепях (позиции IMGT 105-117).
3. Запустить BLAST только для последовательностей с найденным H3.
4. Оставить HSP, которые полностью покрывают CDR-H3 и у которых идентичность H3 попадает в заданный диапазон.

### Выходные файлы в `--output-dir`

- `<seq_hash>.xml` — исходный XML-вывод BLAST.
- `<seq_hash>_filtered_hits.csv` — отфильтрованные попадания для одной последовательности.
- `combined_filtered_hits.csv` — все строки вместе. Столбцы: `source_file, chain_id, hit_id, hit_description, pct_identity, identity, align_length, evalue, bit_score, h3_covered, h3_identity, h3_align_length, h3_pct_identity`.
- `unique_cdr_h3_summary.csv` — столбцы: `cdr_h3_sequence, n_chains_sharing, pdb_ids, seq_hashes, example_sources`.
- `errors.txt` — появляется, если были ошибки.

## normalize_cdr_h3_length.py — только библиотека (без CLI)

Функции и классы для импорта: `LengthBucket`, `max_allowed_edits(length)`, `edit_distance(a, b)`, `is_same_h3(query_h3, hit_h3)`.

Допустимое число правок по длине H3:

| Длина | Допустимо правок |
|---|---|
| 1–4 | 0 |
| 5–8 | 1 |
| 9–14 | 2 |
| 15+ | 3 |

Группа выбирается по более короткой последовательности. Нужен rapidfuzz. Блока `main` для запуска нет.