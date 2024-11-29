#!/bin/bash

cd "$(dirname "$0")"

for python in \
    3.7 \
    3.8 \
    3.9 \
    3.10 \
    3.11 \
    3.12 \
; do
    python3 -m uv pip compile --python-version=$python --allow-unsafe ../setup.py --extra test --extra docs -o "requirements-${python}.txt"
done

wait

ln -sf requirements-3.9.txt requirements.txt

echo "fin"
