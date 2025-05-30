#! /usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import argparse

# parses command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Creates a heatmap from fastANI results")
    parser.add_argument("--input", "-i", help="FastANI results to make a heatmap from")
    parser.add_argument("--out", help="Name of the resulting heatmap")
    return parser.parse_args()

args = parse_args()

# Read data from the file
x = pd.read_table(args.input, header=None)

# Convert to matrix format using pivot
data = x.pivot(index=0, columns=1, values=2)
data.index.name = None
data.columns.name = None

# Convert the percent identity matrix to a distance matrix
# (100% identity => distance 0, 0% identity => distance 1)
distance_matrix = 1 - (data / 100)

# Generate the linkage matrix for the hierarchical clustering
linkage_matrix = linkage(distance_matrix, method='average')

# Plot heatmap with a dendrogram
fig, ax = plt.subplots(figsize=(10, 8))

# Use a clustermap in Seaborn to add a heatmap and dendrogram
clustermap = sns.clustermap(
    data,
    row_linkage=linkage_matrix,
    col_linkage=linkage_matrix,
    cmap='RdYlBu_r',  # Color map for the heatmap
    linewidths=0.5,
    figsize=(20, 16),
    dendrogram_ratio=(.2, .2),  # Adjust the size of the dendrograms
    cbar_pos=(0.02, 0.9, 0.03, 0.15)  # Colorbar position
)

clustermap.fig.suptitle('Heatmap with Cladogram Based on Percent Identity Matrix', y=1.05)

plt.savefig(args.out, bbox_inches="tight")

