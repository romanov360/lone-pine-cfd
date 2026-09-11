#!/usr/bin/env bash
# Reproduce the whole study. Expect a few hours on four cores.
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONPATH=src
mkdir -p results/data results/figures

echo "== verification and validation =="
python3 src/validate.py heated
python3 src/validate.py cavity
python3 src/validate.py cyl

echo "== CFD campaign =="
python3 src/campaign.py grid
python3 src/campaign.py boxsize
python3 src/campaign.py creek
python3 src/campaign.py h2h
python3 src/gci.py || true

echo "== models and sweeps =="
python3 src/props.py
python3 src/geometry.py
python3 src/correlations.py
python3 src/aeration.py
python3 src/axisym_fipy.py
for s in 1 2 3 4 5 6 7 8 9 10 11; do python3 src/sweeps.py "$s"; done

echo "== synthesis =="
python3 src/answer.py | tee results/data/answer.txt

echo "== figures and web export =="
python3 src/figures.py all
python3 src/export_web.py
python3 src/build_report.py
