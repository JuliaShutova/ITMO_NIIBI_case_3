import argparse,json,pandas as pd
p=argparse.ArgumentParser(); p.add_argument("--ranked",required=True); p.add_argument("--metadata",required=True); p.add_argument("--out",required=True); a=p.parse_args()
m=json.load(open(a.metadata)); df=pd.read_csv(a.ranked,sep="\t")
prediction="no_causal_variant_in_vcf"
# Emit a candidate only when the table has a non-empty variant key and score > 0.
# This is a baseline rule; tune threshold on development cases, not hidden cases.
if len(df) and "variant_id" in df and "final_score" in df:
    top=df.iloc[0]
    if pd.notna(top["variant_id"]) and float(top["final_score"])>0:
        v=str(top["variant_id"])
        # Benchmark format excludes assembly prefix.
        prediction=v
with open(a.out,"w") as f:
    f.write("case_id\tprediction\n")
    f.write(f"{m.get('case_id','case')}\t{prediction}\n")
