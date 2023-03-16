#!/usr/local/bin/sh

set -e # exit script on error

# Create and source virtual environment
python3 -m venv .venv
source .venv/bin/activate 

# Install dependencies
pip3 install -r requirements.txt


set +e # return to default shell behaviour 