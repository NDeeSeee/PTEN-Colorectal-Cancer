#!/bin/bash

PYTHON_PATH="/Users/grigoriiandrianov/miniconda3/bin/python"
SCRIPT_DIR="scripts"

 

mkdir -p results

# find linear hotspots
# arguments order: mutations frequency | p-value cutoff
$PYTHON_PATH $SCRIPT_DIR/generate_linear_hotspots.py data/Mutations.tsv 0.005 > results/linear_hotspots.tsv

# find sliding window enrichments
# arguments order: mutations frequency | p-value cutoff
$PYTHON_PATH $SCRIPT_DIR/generate_sliding_windows.py data/Mutations.tsv results/sliding_window.tsv

# find 3D hotspots
# arguments order: 3d interactions | mutations frequency | table with random indicies | linear hotspots | output file
$PYTHON_PATH $SCRIPT_DIR/generate_3d_hotspots.py data/3D_interactions.tsv data/Mutations.tsv data/Random_Table.tsv results/linear_hotspots.tsv results/3dHotspots.tsv

