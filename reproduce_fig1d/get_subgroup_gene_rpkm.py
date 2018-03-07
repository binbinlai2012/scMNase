#!/usr/bin/env python

import os,sys
from string import *

f = open("NIH3T3_rpkm_refGene.f.txt", 'r')
outf1 = open("gene_rpkm10.txt", 'w')
outf2 = open("gene_rpkm1-10.txt", 'w')
outf3 = open("gene_rpkm0.1-1.txt", 'w')
outf4 = open("gene_rpkm0-0.1.txt", 'w')

f.readline()
f.readline()
f.readline()
f.readline()
f.readline()

for line in f:
	tline = line.strip()
	sline = tline.split()
	gn = sline[0]
	rpkm = (float(sline[4]) + float(sline[5]) ) / 2
	
	if rpkm > 10:
		outf1.write("%s\n" % gn)
	elif rpkm > 1:
		outf2.write("%s\n" % gn)
	elif rpkm > 0.1:
		outf3.write("%s\n" % gn)
	else:
		outf4.write("%s\n" % gn)
		
outf1.close()
outf2.close()
outf3.close()
outf4.close()
