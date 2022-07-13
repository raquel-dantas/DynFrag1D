#!/usr/local/bin/sh

set -e # exit script on error

# Create and source virtual environment
python3 -m venv .venv
source .venv/bin/activate 

# Install dependencies
pip3 install -r requirements.txt

# Source akantu built from source 
source ~/projects/lib/akantu/build/akantu_environement.sh


set +e # return to default shell behaviour (You probably want to keep this at the end of this script)