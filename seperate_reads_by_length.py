#!/usr/bin/env python

import os,sys
from string import *

f = open(sys.argv[1], 'r')
outf_80 = open(sys.argv[2], 'w')

for line in f:
	tline = line.strip()
	sline = tline.split()
	l = int(sline[3])
	if l <= 80:
		outf_80.write(line)
f.close()
outf_80.close()

