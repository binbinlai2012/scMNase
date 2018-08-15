#!/usr/bin/env python

import os,sys
from string import *

f = open(sys.argv[1], 'r')
outf_140_180 = open(sys.argv[2], 'w')
#outf_117_139 = open(sys.argv[3], 'w')
#outf_79_116 = open(sys.argv[4], 'w')
#outf_78 = open(sys.argv[5], 'w')
#outf_117 = open(sys.argv[6], 'w')

for line in f:
	tline = line.strip()
	sline = tline.split()
	l = int(sline[3])
	if l >= 140 and l <= 180:
		outf_140_180.write(line)
#	if l >= 117 and l <= 139:
#		outf_117_139.write(line)
#	if l >= 79 and l <= 116:
#		outf_79_116.write(line)
#	if l <= 78:
#		outf_78.write(line)
#	if l >= 117:
#		outf_117.write(line)
f.close()
#outf_140.close()
#outf_117_139.close()
#outf_79_116.close()
#outf_78.close()
#outf_117.close()

