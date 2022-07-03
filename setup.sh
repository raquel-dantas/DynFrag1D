#!/usr/local/bin/sh

python -m venv .venv && \
source .venv/bin/activate && \
pip3 install ipykernel numpy matplotlib progressbar akantu autopep8