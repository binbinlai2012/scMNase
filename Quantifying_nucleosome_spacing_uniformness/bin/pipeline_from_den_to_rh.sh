#!/bin/sh

infile=$1

tag=${infile:0:${#infile}-4}
smfile=${tag}.sm.txt


smooth_density $1 2 $smfile
SBfile=${tag}.sm.SB.txt
GWSBfile=${tag}.sm.gw.SB.txt
GWfile=${tag}.sm.gw.txt

cal_summit_bottom_from_aveprofile -f $smfile -w 0 -o $SBfile
cal_summit_bottom_from_aveprofile -f $smfile -w 1 -o $GWSBfile -g $GWfile

outfile=${tag}.rh.txt
cal_relative_height $GWSBfile $SBfile $outfile
