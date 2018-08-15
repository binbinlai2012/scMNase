#!/usr/bin/env python

import os,sys
from string import *


f = open(sys.argv[1], 'r')

outf = open(sys.argv[2], 'w')
for line in f:
	tline = line.strip()
	sline = tline.split('.')
	outf.write("cal_density_3 %s 500 %s.len_den.txt\n" % (tline, sline[0]) )
outf.close()
f.close()

