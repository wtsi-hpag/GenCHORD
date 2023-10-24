#!/bin/bash

files=( "BC-1-T2T-CCS" "HG002-T2T-CCS" "HG004-T2T-CCS" "OES103-T2T-10X" "OES103-T2T-CCS" "OES117-T2T-CCS" "OES143-T2T-CCS" "OES148-T2T-CCS" "OES152-T2T-CCS" "OES152-T2T-ONT" "HG003-T2T-CCS" "HT-115-T2T-CCS"  "OES103-T2T-blood" "OES103-T2T-HiC" "OES117-T2T-HiC" "OES143-T2T-HiC" "OES148-T2T-HiC" "OES152-T2T-HiC" )
output=( "bc1/bc1-ccs" "hg002/hg002-ccs" "hg004/hg004-ccs" "oes103/oes103-10x" "oes103/oes103-ccs" "oes117/oes117-ccs" "oes143/oes143-ccs" "oes148/oes148-ccs" "oes152/oes152-ccs" "oes152/oes152-ont" "hg003/hg003-ccs" "ht115/ht115-ccs" "oes103/oes103-blood" "oes103/oes103-hic" "oes117/oes117-hic" "oes143/oes143-hic" "oes148/oes148-hic" "oes152/oes152-hic" )
# thinning=( 0.1 0.1 0.1 0.7 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.95 0.95 0.95 0.95 0.95 )
# files=(  "OES148-T2T-HiC" "OES148-T2T-CCS" )
# output=( "oes148-hic" "oes148-ccs" )
# thinning=( 0.99 0.1 )

for thin in 0 0.1 0.5 0.9 0.99 0.999 0.9999; do
	for gamma in 0.1 1 3 5 10; do
		for i in "${!files[@]}"; do
			# echo ${files[i]} " " $thin " " $gamma
			simple="${thin/.}"
			./deforest -f ../Coverage/${files[i]}.dat -sigmaResolution 10 -worker 40 -sigmaMin 5 -sigmaMax 30 -o Output/Thin$simple"Gamma"$gamma/${output[i]} -L 1000000 -gamma $gamma -thin 10 -accelerate 50 -Qmax 14 -alpha 0.01  -smooth $thin
		done
	done
done