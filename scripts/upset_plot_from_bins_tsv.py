#! /usr/bin/env python3

import argparse
from upsetplot import UpSet
from upsetplot import from_contents
import matplotlib.pyplot as plt
import matplotlib
from csv import DictReader
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

"""
Orthogroup	Kurthia_huakuii	Solibacillus_kalamii	Kurthia_zopfii	Kurthia_sibirica	PROKKA_10142021	Total
OG0000000	1	0	5	9	19	34
OG0000001	22	0	1	2	2	27
OG0000002	5	5	4	7	6	27
OG0000003	3	20	1	1	1	26
"""

# parses command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Creates an Upset plot from orthogroup data")
    parser.add_argument("--genes", "-g", help="Orthogroup gene count csv file")
    parser.add_argument("--out", default="orthogroup_upsetplot.png", help="Name of Upset plot file")
    return parser.parse_args()

args = parse_args()

contents_dict = {}
# import binning tsv into a contents dictionary
bin_dict = DictReader(open(args.genes), delimiter="\t")

for row in bin_dict:
    curr_bin = ""
    for sample, count in row.items():
        # if its the Total row, ignore as sample
        if sample == "Total":
            continue
        # if sample is orthogroup, count is actually the sequence of interest
        if sample == "Orthogroup":
            contents_dict.setdefault(count, [])
            curr_bin = count
            continue
        # otherwise, sample is valid and will have count to add to contents
        if float(count) > 0:
            contents_dict[curr_bin].append(sample)

# transpose dictionary contents
transpose_contents = {}
for og, sample_list in contents_dict.items():
    for sample in sample_list:
        transpose_contents.setdefault(sample, [])
        transpose_contents[sample].append(og)

orthogroups = from_contents(transpose_contents)
matplotlib.rcParams["font.size"] = 7.25
ax_dict = UpSet(orthogroups, subset_size="count", show_counts="%d", show_percentages=True, sort_categories_by="-input", sort_by="cardinality").plot()

# Save the plot to a file
plt.savefig(args.out)

# Show the plot
#plt.show()
