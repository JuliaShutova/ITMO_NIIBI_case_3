import json,sys
src,out=sys.argv[1:3]
with open(src,encoding="utf-8") as f: d=json.load(f)
required=["case_id","genome_assembly"]
missing=[x for x in required if not d.get(x)]
if missing: raise SystemExit("Missing required metadata: "+", ".join(missing))
with open(out,"w") as f: json.dump(d,f,indent=2)
