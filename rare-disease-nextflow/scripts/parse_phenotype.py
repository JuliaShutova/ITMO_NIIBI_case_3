import argparse,json
p=argparse.ArgumentParser()
p.add_argument("input"); p.add_argument("--parsed",required=True); p.add_argument("--present",required=True)
p.add_argument("--excluded",required=True); p.add_argument("--phenopacket",required=True)
a=p.parse_args()
with open(a.input,encoding="utf-8") as f: d=json.load(f)
present=[]; excluded=[]; features=[]
for item in d.get("phenotypicFeatures",[]):
    term=item.get("type",{}); hid=term.get("id"); label=term.get("label")
    if not hid: continue
    neg=bool(item.get("excluded",False))
    (excluded if neg else present).append(hid)
    features.append({"id":hid,"label":label,"excluded":neg})
parsed={"case_id":d.get("id"),"sex":d.get("sex","UNKNOWN_SEX"),
        "present_hpo":sorted(set(present)),"excluded_hpo":sorted(set(excluded)),
        "phenotype_truncated":d.get("phenotype_truncated",False),"features":features}
for path,vals in [(a.present,present),(a.excluded,excluded)]:
    with open(path,"w") as f: f.write("\n".join(sorted(set(vals)))+("\n" if vals else ""))
with open(a.parsed,"w") as f: json.dump(parsed,f,indent=2)
# Phenopacket-like JSON for LIRICAL. The actual accepted schema is version-specific;
# validate this file against the installed LIRICAL release before production use.
packet={"id":d.get("id","case"),"subject":{"id":d.get("id","case")},
        "phenotypicFeatures":[{"type":{"id":x},"excluded":False} for x in sorted(set(present))]+
                             [{"type":{"id":x},"excluded":True} for x in sorted(set(excluded))]}
with open(a.phenopacket,"w") as f: json.dump(packet,f,indent=2)
