import argparse, gzip, csv
p=argparse.ArgumentParser(); p.add_argument('--vcf',required=True); p.add_argument('--out',required=True); a=p.parse_args()
op=gzip.open if a.vcf.endswith('.gz') else open
rows=[]
with op(a.vcf,'rt') as f:
    sample_idx=None
    for line in f:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            h=line.rstrip().split('\t'); sample_idx=9 if len(h)>9 else None; continue
        if line.startswith('#') or not line.strip(): continue
        c=line.rstrip().split('\t')
        gt='NA'
        if sample_idx is not None and len(c)>sample_idx:
            gt=c[sample_idx].split(':')[0]
        rows.append({'chrom':c[0],'pos':c[1],'ref':c[3],'alt':c[4], 'genotype':gt,
                     'variant_id':f'{c[0]}:{c[1]}:{c[3]}:{c[4]}',
                     'variant_type':('SNV' if len(c[3])==1 and len(c[4])==1 else ('InDel' if len(c[3])<=50 and len(c[4])<=50 else 'SV')),
                     'gene':'','gene_id':'','consequence':'','impact':'','hgvs_c':'','hgvs_p':'',
                     'clinvar':'NA','gnomad_af':'NA'})
fields=['chrom','pos','ref','alt','genotype','variant_id','variant_type','gene','gene_id','consequence','impact','hgvs_c','hgvs_p','clinvar','gnomad_af']
with open(a.out,'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)
