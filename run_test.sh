#!/bin/bash
# Install toolbox package locally in the job's working directory
# export PYTHONPATH=$PWD:$PYTHONPATH
pip install --user -e .
# create test folder if it doesn't exist
mkdir -p test
# Run your script
python3 multiple_patterns.py init/7_0/ test/