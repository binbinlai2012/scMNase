#!/usr/bin/env python

import os,sys
from string import *
import math, math
import scipy, scipy.stats
from scipy.stats import hypergeom

backgroundfile = sys.argv[1]
obsfile = sys.argv[2]

outf1 = open(sys.argv[3], 'w')
outf2 = open(sys.argv[4], 'w')

def readr( infile ):
	f = open(infile, 'r')
	cell_r = {}
	cell_count = {}
	for line in f:
		tline = line.strip()
		sline = tline.split()
		c = sline[0]
		r = float(sline[1])
		n = int(sline[2])
		tn = int(sline[3])
		cell_r.update({c:r})
		cell_count.update({c:(n,tn)})
	f.close()
	return (cell_r, cell_count )

(cell_r_bg, cell_count_bg) = readr(backgroundfile)
(cell_r_obs, cell_count_obs) = readr(obsfile)

cells = []
scs = []
for cell in cell_r_obs.keys():
	if cell in cell_r_bg.keys():
		r_obs = cell_r_obs[cell]
		r_bg = cell_r_bg[cell]
		(n_obs, tn_obs) = cell_count_obs[cell]
		(n_bg, tn_bg) = cell_count_bg[cell]
		
	#	print cell, r_obs, r_bg
		sc = math.log(r_obs/r_bg, 2)
	#	sc = r_bg - r_obs
	#	print sc
		p = hypergeom.cdf(n_obs, tn_bg, n_bg, tn_obs, loc=0)
		cells.append(cell)
		scs.append(sc)
		outf1.write( "%s\t%f\t%f\n" % (cell, sc, p))
outf1.close()

zscs = scipy.stats.zscore(scs)
print len(scs), len(cells)
for i in range(len(cells)):
	outf2.write("%s\t%f\n" % (cells[i], zscs[i]))
	
outf2.close()

		
