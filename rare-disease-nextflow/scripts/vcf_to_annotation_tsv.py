import argparse,gzip,csv
p=argparse.ArgumentParser(); p.add_argument("--vcf",required=True); p.add_argument("--out",required=True); a=p.parse_args()
op=gzip.open if a.vcf.endswith(".gz") else open
rows=[]
with op(a.vcf,"rt") as f:
    for line in f:
        if line.startswith("#"): continue
        c=line.rstrip().split("\t"); info={}
        for x in c[7].split(";"):
            if "=" in x:
                k,v=x.split("=",1); info[k]=v
        anns=info.get("ANN","").split(",") if info.get("ANN") else [""]
        for ann in anns:
            z=ann.split("|")
            rows.append({"chrom":c[0],"pos":c[1],"ref":c[3],"alt":c[4],
                         "genotype":c[9].split(":")[0] if len(c)>9 else "NA",
                         "variant_id":f"{c[0]}:{c[1]}:{c[3]}:{c[4]}",
                         "ann_allele":z[0] if len(z)>0 else "",
                         "consequence":z[1] if len(z)>1 else "",
                         "impact":z[2] if len(z)>2 else "",
                         "gene":z[3] if len(z)>3 else "",
                         "gene_id":z[4] if len(z)>4 else "",
                         "feature":z[6] if len(z)>6 else "",
                         "hgvs_c":z[9] if len(z)>9 else "",
                         "hgvs_p":z[10] if len(z)>10 else "",
                         "clinvar":info.get("CLNSIG","NA"),
                         "gnomad_af":info.get("gnomAD_AF","NA")})
fields=list(rows[0]) if rows else ["chrom","pos","ref","alt","genotype","variant_id","gene","gene_id","consequence","impact","hgvs_c","hgvs_p","clinvar","gnomad_af"]
with open(a.out,"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerows(rows)
