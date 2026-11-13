import json, tempfile, subprocess, pathlib, unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
class PipelineHelpers(unittest.TestCase):
    def test_parse_phenotype_keeps_excluded_separate(self):
        with tempfile.TemporaryDirectory() as d:
            d=pathlib.Path(d)
            src=d/"p.json"; src.write_text(json.dumps({"id":"x","phenotypicFeatures":[
                {"type":{"id":"HP:0000001","label":"A"}},
                {"type":{"id":"HP:0000002","label":"B"},"excluded":True}]}))
            subprocess.run(["python",str(ROOT/"scripts/parse_phenotype.py"),str(src),
                "--parsed",str(d/"o.json"),"--present",str(d/"present"),
                "--excluded",str(d/"excluded"),"--phenopacket",str(d/"packet")],check=True)
            self.assertEqual((d/"present").read_text().strip(),"HP:0000001")
            self.assertEqual((d/"excluded").read_text().strip(),"HP:0000002")
if __name__=="__main__": unittest.main()
