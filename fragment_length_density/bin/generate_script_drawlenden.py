#!/usr/bin/env python

import os,sys
from string import *


f = open(sys.argv[1], 'r')

outf = open(sys.argv[2], 'w')
for line in f:
	tline = line.strip()
	sline = tline.split('.')
	outf.write("Rscript drawprofile.r %s.len_den.txt %s.len_den.png\n" % (sline[0], sline[0]) )
outf.close()
f.close()

