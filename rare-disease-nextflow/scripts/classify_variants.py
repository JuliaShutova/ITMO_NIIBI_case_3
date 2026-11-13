import argparse,gzip,re,subprocess
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument("--vcf",required=True); p.add_argument("--outdir",required=True); p.add_argument("--report",required=True)
a=p.parse_args(); out=Path(a.outdir); out.mkdir(exist_ok=True,parents=True)
groups={"snv":[],"indel":[],"sv":[],"unclassified":[]}
header=[]
with gzip.open(a.vcf,"rt") if a.vcf.endswith(".gz") else open(a.vcf) as f:
    for line in f:
        if line.startswith("#"): header.append(line); continue
        c=line.rstrip("\n").split("\t")
        if len(c)<8: continue
        ref,alt=c[3],c[4]
        if "," in alt: # split multi-allelic records upstream with bcftools norm -m -any
            groups["unclassified"].append(line); continue
        if alt.startswith("<") or "[" in alt or "]" in alt or alt=="*":
            groups["sv"].append(line); continue
        if len(ref)==1 and len(alt)==1: groups["snv"].append(line)
        elif abs(len(ref)-len(alt))<=50: groups["indel"].append(line)
        elif abs(len(ref)-len(alt))>50: groups["sv"].append(line)
        else: groups["unclassified"].append(line)
with open(a.report,"w") as rep:
    rep.write("category\tcount\tfile\n")
    for cat,rows in groups.items():
        rep.write(f"{cat}\t{len(rows)}\t{cat}.vcf.gz\n")
        if rows:
            tmp=out/f"{cat}.vcf"
            tmp.write_text("".join(header+rows))
            subprocess.run(["bgzip","-f",str(tmp)],check=True)
            subprocess.run(["tabix","-f","-p","vcf",str(tmp)+".gz"],check=True)
