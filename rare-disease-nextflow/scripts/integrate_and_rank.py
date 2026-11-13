import argparse,json,glob,os
from pathlib import Path
import pandas as pd
p=argparse.ArgumentParser()
for x in ["small","sv","phen2gene","lirical","base","classification","metadata","phenotype","out","integrated","report"]:
    p.add_argument("--"+x,required=True)
a=p.parse_args()
def read_many(spec):
    paths=[]
    for token in spec.split():
        paths.extend(glob.glob(token))
    frames=[]
    for x in paths:
        if os.path.isfile(x) and os.path.getsize(x)>0:
            try: frames.append(pd.read_csv(x,sep="\t",comment="!",dtype=str))
            except Exception: pass
    return pd.concat(frames,ignore_index=True,sort=False) if frames else pd.DataFrame()
small=read_many(a.small); sv=read_many(a.sv); base=read_many(a.base); variants=pd.concat([small,sv,base],ignore_index=True,sort=False)
if not variants.empty and "variant_id" in variants.columns:
    variants=variants.drop_duplicates(subset=["variant_id"],keep="first")
p2g=read_many(a.phen2gene); lir=read_many(a.lirical)
# Harmonize gene symbols where possible; keep raw phenotype outputs available separately.
if not variants.empty:
    variants["gene"]=variants.get("gene",pd.Series("",index=variants.index)).fillna("")
    if not p2g.empty:
        gene_col=next((c for c in p2g.columns if c.lower() in {"gene","gene_symbol","symbol","genesymbol"}),None)
        score_col=next((c for c in p2g.columns if "score" in c.lower()),None)
        rank_col=next((c for c in p2g.columns if "rank" in c.lower()),None)
        if gene_col:
            p2g2=p2g[[gene_col]+([score_col] if score_col else [])+([rank_col] if rank_col else [])].copy()
            p2g2=p2g2.rename(columns={gene_col:"gene", **({score_col:"phen2gene_score"} if score_col else {}), **({rank_col:"phen2gene_rank"} if rank_col else {})})
            variants=variants.merge(p2g2.drop_duplicates("gene"),on="gene",how="left")
    # Conservative transparent heuristic. Do not claim calibrated pathogenicity probability.
    consequence=variants.get("consequence",pd.Series("",index=variants.index)).fillna("").str.lower()
    impact=variants.get("impact",pd.Series("",index=variants.index)).fillna("").str.upper()
    variants["functional_score"]=impact.map({"HIGH":3,"MODERATE":2,"LOW":1,"MODIFIER":0}).fillna(0)
    variants.loc[consequence.str.contains("stop_gained|frameshift|splice_acceptor|splice_donor",regex=True),"functional_score"]=4
    clin=variants.get("clinvar",pd.Series("",index=variants.index)).fillna("").str.lower()
    variants["clinical_score"]=clin.str.contains("pathogenic").astype(int)*4
    af=pd.to_numeric(variants.get("gnomad_af",pd.Series(index=variants.index,dtype=float)),errors="coerce")
    variants["population_score"]=0
    variants.loc[af.notna() & (af<0.0001),"population_score"]=2
    variants.loc[af.notna() & (af>=0.01),"population_score"]=-2
    ps=pd.to_numeric(variants.get("phen2gene_score",pd.Series(index=variants.index,dtype=float)),errors="coerce")
    variants["phenotype_score"]=ps.fillna(0)
    variants["final_score"]=variants[["functional_score","clinical_score","population_score","phenotype_score"]].sum(axis=1)
    variants["evidence_for"]=variants.apply(lambda r:";".join([x for x in [
        "high/moderate functional impact" if r.get("functional_score",0)>=2 else "",
        "ClinVar pathogenic assertion" if r.get("clinical_score",0)>0 else "",
        "rare in gnomAD (if observed)" if r.get("population_score",0)>0 else "",
        "Phen2Gene support" if pd.notna(r.get("phen2gene_score",None)) else ""
    ] if x]),axis=1)
    variants["evidence_against"]=variants.apply(lambda r:"common population frequency" if r.get("population_score",0)<0 else "",axis=1)
    variants["explanation"]=variants.apply(lambda r:f"Evidence-based heuristic score={r.get('final_score',0)}; not a calibrated probability.",axis=1)
    variants=variants.sort_values(["final_score"],ascending=False).reset_index(drop=True)
    variants.insert(0,"rank",range(1,len(variants)+1))
else:
    variants=pd.DataFrame(columns=["rank","variant_id","chrom","pos","ref","alt","gene","final_score","explanation"])
variants.to_csv(a.out,sep="\t",index=False)
variants.to_csv(a.integrated,sep="\t",index=False)
meta=json.load(open(a.metadata)); phenotype=json.load(open(a.phenotype))
report={"case_id":meta.get("case_id"),"assembly":meta.get("genome_assembly"),
        "n_ranked_variants":int(len(variants)),"n_phen2gene_rows":int(len(p2g)),
        "n_lirical_rows":int(len(lir)),
        "warning":"Baseline heuristic only. Missing annotations are not benign evidence. Disease-level LIRICAL scores are not automatically variant-level scores."}
json.dump(report,open(a.report,"w"),indent=2)
