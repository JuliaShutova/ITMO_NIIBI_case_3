import argparse, json, csv, subprocess, sys
from pathlib import Path

p=argparse.ArgumentParser()
for x in ("metadata","phenotype","variants"): p.add_argument("--"+x, required=True)
p.add_argument("--out", required=True)
a=p.parse_args()
report={"status":"PASS","errors":[],"warnings":[],"files":{}}
def load_json(path):
    with open(path, encoding="utf-8") as f: return json.load(f)
for name in ("metadata","phenotype","variants"):
    path=Path(getattr(a,name))
    report["files"][name]={"exists":path.exists(),"readable":path.is_file()}
    if not path.is_file(): report["errors"].append(f"{name}: file missing"); continue
    try:
        if name in ("metadata","phenotype"): load_json(path)
        else:
            with path.open(encoding="utf-8-sig") as f: first=f.readline().strip()
            if first.startswith("##fileformat=VCF"):
                r=subprocess.run(["bcftools","view","-h",str(path)],capture_output=True,text=True)
                if r.returncode: report["errors"].append("bcftools cannot read VCF: "+r.stderr[-1000:])
            elif "\t" in first and {"chrom","pos","ref","alt"}.issubset(set(first.lower().split("\t"))):
                report["warnings"].append("variants input detected as tabular TSV, not standard VCF")
            else: report["errors"].append("variants file is neither recognizable VCF nor expected TSV")
    except Exception as e: report["errors"].append(f"{name}: {e}")
try:
    m=load_json(a.metadata); h=load_json(a.phenotype)
    if m.get("case_id") and h.get("id") and m["case_id"] != h["id"]:
        report["errors"].append("case_id in metadata does not match phenotype.id")
    if m.get("genome_assembly") != "GRCh38":
        report["warnings"].append(f"assembly is {m.get('genome_assembly')}; verify all databases match")
except Exception: pass
if report["errors"]: report["status"]="FAIL"
with open(a.out,"w") as f: json.dump(report,f,indent=2)
if report["errors"]: sys.exit(1)
