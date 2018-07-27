#!/bin/sh

infile=$1

spar=$2


tag=${infile:0:${#infile}-4}
smfile=${tag}.sm.txt


#smooth_density $1 2 $smfile
splinefile=${tag}.spline.txt
Rscript ${PWD}/spline.r $infile $spar $splinefile

SBfile=${tag}.spline.SB.txt
GWSBfile=$SBfile


cal_summit_bottom_from_aveprofile -f $splinefile -w 0 -o $SBfile


outfile=${tag}.rh.txt
cal_relative_height $GWSBfile $SBfile $outfile

