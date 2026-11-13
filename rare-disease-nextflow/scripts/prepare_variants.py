import argparse,csv,json,re,sys
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument("--input",required=True); p.add_argument("--metadata",required=True); p.add_argument("--out",required=True)
a=p.parse_args()
meta=json.load(open(a.metadata)); assembly=meta.get("genome_assembly","GRCh38")
src=Path(a.input)
with src.open(encoding="utf-8-sig") as f: first=f.readline().strip()
if first.startswith("##fileformat=VCF"):
    # Preserve original records and header. Add a sample header only if already valid;
    # bcftools performs strict validation in the next process.
    Path(a.out).write_text(src.read_text())
else:
    with src.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        required={"chrom","pos","ref","alt","genotype"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise SystemExit(f"TSV must contain columns: {sorted(required)}")
        rows=list(reader)
    sample=meta.get("case_id","SAMPLE")
    with open(a.out,"w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##reference={assembly}\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+sample+"\n")
        for i,r in enumerate(rows,1):
            chrom=r["chrom"]; chrom=chrom if chrom.startswith("chr") else "chr"+chrom
            pos=int(r["pos"]); ref=r["ref"]; alt=r["alt"]; gt=r["genotype"]
            if pos<1 or not ref or not alt or gt not in {"0/0","0/1","1/0","1/1","./.","0|0","0|1","1|0","1|1"}:
                raise SystemExit(f"Invalid variant row {i}: {r}")
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gt}\n")
