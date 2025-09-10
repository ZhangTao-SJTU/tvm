#!/bin/bash
# Install toolbox package locally in the job's working directory
#export PYTHONPATH=$PWD:$PYTHONPATH
source ~/.bashrc
pip install --user -e .
# create main run folder if it doesn't exist
mkdir -p multiple_patterns_2_cells
# Run your script
python multiple_patterns.py init/7_1/ multiple_patterns_2_cells/test/
