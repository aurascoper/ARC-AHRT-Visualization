#!/bin/bash
# Download the Geant4 11.3.2 dataset set (manifest-extracted), verify md5, extract.
set -e
D=~/Developer/AHRT/topas/G4Data-11.3.2
mkdir -p "$D"; cd "$D"
python3 - <<'EOF' > /tmp/g4data_dl.sh
import json
rows = json.load(open("/Users/aurascoper/Developer/AHRT/topas/g4data_1132_manifest.json"))
for r in rows:
    fn = r["FILENAME"]; ver = r["VERSION"]; ext = r["EXTENSION"]
    url = f"https://cern.ch/geant4-data/datasets/{fn}.{ver}.{ext}"
    print(f'curl -sfL --retry 3 -o "{fn}.{ver}.{ext}" "{url}" && md5 -q "{fn}.{ver}.{ext}" | grep -qi {r["MD5SUM"]} && tar xzf "{fn}.{ver}.{ext}" && echo "OK {fn} {ver}" || echo "FAIL {fn} {ver}"')
EOF
bash /tmp/g4data_dl.sh