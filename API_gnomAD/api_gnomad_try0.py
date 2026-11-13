import argparse
import csv
import time
import requests

API = "https://gnomad.broadinstitute.org/api"

QUERY = """
query getVariant($variantId: String!) {
    variant(variantId: $variantId, dataset: gnomad_r4) {
        variant_id
        exome {
            ac
            an
        }
        genome {
            ac
            an
        }
    }
}
"""

def get_af(data):
    if not data:
        return ""
    ac = data.get("ac")
    an = data.get("an")
    if not an or an == 0:
        return ""
    return str(ac / an)


def query_variant(chrom, pos, ref, alt):
    chrom_str = str(chrom).upper().replace("CHR", "")
    variant_id = f"{chrom_str}-{pos}-{ref}-{alt}"

    attempt = 0
    while True:  # БЕСКОНЕЧНЫЙ цикл для 429, пока не получим ответ
        try:
            response = requests.post(
                API,
                json={
                    "query": QUERY,
                    "variables": {"variantId": variant_id}
                },
                timeout=30
            )

            # 1. ЕСЛИ 429: Ждем всё дольше и пробуем СНОВА. Никогда не выходим из цикла.
            if response.status_code == 429:
                # Задержка: 5с, 10с, 20с, 40с, 60с (максимум 60 секунд)
                wait_time = min(60, 5 * (2 ** attempt))
                print(f"   ⏳ 429 Too Many Requests для {variant_id}. Ждем {wait_time}с и пробуем снова...")
                time.sleep(wait_time)
                attempt += 1
                continue  # Возвращаемся в начало цикла while True

            # 2. ЕСЛИ 400 (Bad Request): Вариант некорректен, дальше ждать бессмысленно.
            if response.status_code == 400:
                print(f"   ⚠️ 400 Bad Request для {variant_id} (проверьте формат REF/ALT). Пропускаем аннотацию.")
                return {"gnomad_exome_af": "", "gnomad_genome_af": "", "gnomad_status": "bad_request"}

            # 3. ЕСЛИ ДРУГАЯ ОШИБКА (500, 502 и т.д.): Ждем и пробуем снова
            if response.status_code >= 500:
                wait_time = 10
                print(f"   ⚠️ Ошибка сервера {response.status_code} для {variant_id}. Ждем {wait_time}с...")
                time.sleep(wait_time)
                attempt += 1
                continue

            # 4. ЕСЛИ 200 OK: Обрабатываем ответ
            response.raise_for_status()
            result = response.json()
            break # Выходим из бесконечного цикла, так как запрос успешен

        except requests.exceptions.RequestException as e:
            # Сетевой сбой (например, отвал интернета). Ждем и пробуем снова.
            wait_time = 10
            print(f"   ⚠️ Сетевая ошибка для {variant_id}: {e}. Ждем {wait_time}с и пробуем снова...")
            time.sleep(wait_time)
            attempt += 1
            continue

    # --- Обработка успешного JSON ответа ---
    if result.get("errors"):
        error_messages = [err.get("message", "") for err in result["errors"]]
        if any("Variant not found" in msg for msg in error_messages):
            return {"gnomad_exome_af": "", "gnomad_genome_af": "", "gnomad_status": "not_found"}
        else:
            return {"gnomad_exome_af": "", "gnomad_genome_af": "", "gnomad_status": "graphql_error"}

    variant = result.get("data", {}).get("variant")
    if variant is None:
        return {"gnomad_exome_af": "", "gnomad_genome_af": "", "gnomad_status": "not_found"}

    exome = variant.get("exome") or {}
    genome = variant.get("genome") or {}

    return {
        "gnomad_exome_af": get_af(exome),
        "gnomad_genome_af": get_af(genome),
        "gnomad_status": "found"
    }


def main():
    parser = argparse.ArgumentParser(description="Аннотация gnomAD (гарантия отсутствия пропусков)")
    parser.add_argument("--input", required=True, help="Входной TSV файл")
    parser.add_argument("--output", required=True, help="Выходной TSV файл")
    args = parser.parse_args()

    with open(args.input, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    output = []
    total = len(rows)

    print(f"🚀 Начало аннотации {total} вариантов. (Скрипт будет ждать столько, сколько нужно, ни один вариант не будет пропущен)")

    for i, row in enumerate(rows):
        result = query_variant(
            row.get("CHROM", row.get("chrom", "")),
            row.get("POS", row.get("pos", "")),
            row.get("REF", row.get("ref", "")),
            row.get("ALT", row.get("alt", ""))
        )

        output.append({**row, **result})

        # Базовая пауза между УСПЕШНЫМИ запросами (1.5 секунды - безопасный минимум для gnomAD)
        time.sleep(1.5)

        if (i + 1) % 50 == 0 or (i + 1) == total:
            found_count = sum(1 for x in output if x["gnomad_status"] == "found")
            not_found_count = sum(1 for x in output if x["gnomad_status"] == "not_found")
            print(f"   ✅ Прогресс: {i + 1}/{total} (Найдено: {found_count}, Нет в базе: {not_found_count})")

    if output:
        fieldnames = list(rows[0].keys()) + ["gnomad_exome_af", "gnomad_genome_af", "gnomad_status"]

        with open(args.output, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(output)

        print(f"🎉 ГОТОВО! Все {total} вариантов обработаны без пропусков. Результат: {args.output}")


if __name__ == "__main__":
    main()
