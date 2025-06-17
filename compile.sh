#!/bin/bash

# This script empties the build/ folder and compiles the project

# To run this script from a terminal:
# 1. Navigate to the project's root folder (where this script should be located...)
# 2. Run the following command in the terminal: chmod +x compile.sh && ./compile.sh

# Create the build/ folder if it doesn't exist
if [ ! -d "build" ]; then
    mkdir build
    fi
# Empty the build/ folder
rm -rf build/*
# Compile the project
cd build/ && cmake ../ && make