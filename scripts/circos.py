#! /usr/bin/env python

from pycirclize import Circos
from pycirclize.parser import Gff
from pycirclize.parser import Genbank
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from Bio.SeqUtils import gc_fraction
import argparse

# parses command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Creates a Circos plot from an assembly")
    parser.add_argument("--species", "-s", help="isolate species")
    parser.add_argument("--gff", help="gff file from PROKKA")
    parser.add_argument("--gbk", help="genbank file from PROKKA")
    parser.add_argument("--bed", help="bed file with per base coverage")
    return parser.parse_args()

args = parse_args()

# Load GFF and GBK files and parse
gff_file = args.gff
gbk_file = args.gbk
gff = Gff(gff_file)
gbk = Genbank(gbk_file)

# Load BED file and parse
bed_file = args.bed
cov = {}
total_cov = 0
with open(bed_file) as bed:
	for line in bed:
		contig, pos, coverage = line.rstrip().split("\t")
		total_cov += int(coverage)
		cov.setdefault(contig, {})
		cov[contig][int(pos)] = int(coverage)
cov_list = list(cov.values())

# track which contig is at n50
half_genome_length = gff.full_genome_length/2
total_length = 0
n50_marked = False

# get whole genome gc content
gc = round(100 * gc_fraction(gbk.full_genome_seq, ambiguous="ignore"), 2)

# Get contig genome seqid & size, features dict
seqid2size = gff.get_seqid2size()
seqid2features = gff.get_seqid2features(feature_type=None)
circos = Circos(seqid2size, space=1.25)

for idx, sector in enumerate(circos.sectors):
    # Plot outer track
	outer_track = sector.add_track((98, 100))
	outer_track.axis(fc="white")
	major_interval = 100000
	minor_interval = int(major_interval / 10)
	if sector.size > minor_interval:
		outer_track.xticks_by_interval(major_interval, label_formatter=lambda v: f"{v / 1000:.0f} Kb", 
										label_orientation="vertical", 
										show_label=True if idx==0 else False)
		outer_track.xticks_by_interval(minor_interval, tick_length=1, show_label=False)

	# Plot forward/reverse CDS, rRNA, tRNA tracks
	f_cds_track = sector.add_track((91, 98), r_pad_ratio=0.1)
	r_cds_track = sector.add_track((84, 91), r_pad_ratio=0.1)
	rna_track = sector.add_track((77, 84), r_pad_ratio=0.1)
	for feature in seqid2features[sector.name]:
		if feature.type == "CDS":
			if feature.location.strand == 1:
				f_cds_track.genomic_features([feature], fc="tomato")
			else:
				r_cds_track.genomic_features([feature], fc="skyblue")
		elif feature.type == "rRNA":
			rna_track.genomic_features([feature], fc="forestgreen")
		elif feature.type == "tRNA":
			rna_track.genomic_features([feature], fc="purple")
		
	# for marking n50
	total_length += sector.size
	if (total_length >= half_genome_length) and not n50_marked:
		n50 = sector.size
		sector.text("N50")
		n50_marked = True

	# Plot GC content
	try:
		gc_content_track = sector.add_track((62, 77), r_pad_ratio=0.1)

		pos_list, raw_gc_content = gbk.calc_gc_content(seq=gbk.records[idx].seq)
		gc_contents = raw_gc_content - gc
		positive_gc_contents = np.where(gc_contents > 0, gc_contents, 0)
		negative_gc_contents = np.where(gc_contents < 0, gc_contents, 0)
		abs_max_gc_content = np.max(np.abs(gc_contents))
		vmin, vmax = -abs_max_gc_content, abs_max_gc_content

		gc_content_track.fill_between(
			pos_list, positive_gc_contents, 0, vmin=vmin, vmax=vmax, color="black"
		)
		gc_content_track.fill_between(
			pos_list, negative_gc_contents, 0, vmin=vmin, vmax=vmax, color="grey"
		)
	except ValueError:
		pass

	# plot coverage
	try:
		cov_track = sector.add_track((47, 62), r_pad_ratio=0.1)
		pos_list, coverages = zip(*cov_list[idx].items())
		cov_track.fill_between(pos_list, coverages, color="orange")
	except ValueError:
		pass

# get species from user as input if not in command line
if args.species is None:
	species = input("Enter your isolate species: ")
else:
	species = args.species

circos.text(f"{species}\nAssembly Length: {gff.full_genome_length:,}\nContigs: {len(circos.sectors):,}\nGC Content: {gc}%\nN50: {n50:,}\nCoverage: {round(total_cov/gff.full_genome_length, 1):,}", r=18, size=8)

fig = circos.plotfig()
_ = circos.ax.legend(
	handles=[
		Patch(color="tomato", label="Forward CDS"),
		Patch(color="skyblue", label="Reverse CDS"),
		Patch(color="forestgreen", label="rRNA"),
		Patch(color="purple", label="tRNA"),
		Patch(color="orange", label="Coverage"),
		Line2D([], [], color="black", label="Positive GC Content", marker="^", ms=6, ls="None"),
    	Line2D([], [], color="grey", label="Negative GC Content", marker="v", ms=6, ls="None")
	],
	bbox_to_anchor=(0.5, 0.4),
	loc="center",
	ncols=2,
	fontsize=8
)

fig.savefig("circos_plot.png", dpi=600)
#plt.show()
