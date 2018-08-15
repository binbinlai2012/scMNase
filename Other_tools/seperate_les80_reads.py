#!/usr/bin/env python

import os,sys
from string import *

f = open(sys.argv[1], 'r')
#outf_100_180 = open(sys.argv[2], 'w')
#outf_117_139 = open(sys.argv[3], 'w')
#outf_79_116 = open(sys.argv[4], 'w')
outf_80 = open(sys.argv[2], 'w')
#outf_117 = open(sys.argv[6], 'w')

for line in f:
	tline = line.strip()
	sline = tline.split()
	l = int(sline[3])
#	if l >= 100 and l <= 180:
#		outf_100_180.write(line)
#	if l >= 117 and l <= 139:
#		outf_117_139.write(line)
#	if l >= 79 and l <= 116:
#		outf_79_116.write(line)
	if l <= 80:
		outf_80.write(line)
#	if l >= 117:
#		outf_117.write(line)
f.close()
#outf_140.close()
#outf_117_139.close()
#outf_79_116.close()
outf_80.close()
#outf_117.close()

