#!/usr/bin/env python

import os,sys
from string import *

f = open(sys.argv[1], 'r')
outf1 = open(sys.argv[2], 'w')
outf2 = open(sys.argv[3], 'w')
prim_c = 0
noprim_c = 0
for line in f:
	tline = line.strip()
	sline = tline.split()
	cell = sline[0]
	sc = float(sline[1])
	p = float(sline[2])
	if p <= 0.05 and sc < -0.1:
		outf1.write("%s\n" % (cell))
		prim_c += 1
	else:
		outf2.write("%s\n" % (cell))
		noprim_c += 1
outf1.close()
outf2.close()

f.close()
print prim_c, noprim_c
