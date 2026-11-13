import argparse,glob,pandas as pd,sys
p=argparse.ArgumentParser(); p.add_argument("--input-glob",required=True); p.add_argument("--out",required=True); a=p.parse_args()
matches=glob.glob(a.input_glob)
if not matches: raise SystemExit("Phen2Gene output not found; check the installed CLI/version and output prefix")
# Preserve tool output as tabular data; normalize only common rank/score/gene columns.
f=matches[0]
try: df=pd.read_csv(f,sep=None,engine="python",comment="#")
except Exception as e: raise SystemExit(f"Cannot parse Phen2Gene output {f}: {e}")
df.to_csv(a.out,sep="\t",index=False)
