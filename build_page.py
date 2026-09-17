#!/usr/bin/env python3
"""Assemble web/template.html + data/build.json + web/coast.min.json -> dist/index.html"""
import json, os, sys
root = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(root, "data/build.json")))
coast = open(os.path.join(root, "web/coast.min.json")).read()
block = (f"const META={json.dumps(d['meta'])};\n"
         f'const B64="{d["ocean_b64"]}";\n'
         f"const COAST={coast};\n"
         f"const ENSO={json.dumps(d['enso'])};\n"
         f'const LB64="{d["land_b64"]}";\n'
         f"const SURV={json.dumps(d.get('surv') or {})};\n")
tpl = open(os.path.join(root, "web/template.html")).read()
assert "/*__DATA__*/" in tpl
out = tpl.replace("/*__DATA__*/", block)
os.makedirs(os.path.join(root, "dist"), exist_ok=True)
open(os.path.join(root, "dist/index.html"), "w").write(out)
print("dist/index.html", round(len(out)/1e6, 2), "Mo")
