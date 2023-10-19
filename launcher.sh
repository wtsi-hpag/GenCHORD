#!/bin/bash

# files=( "BC-1-T2T-CCS" "HG002-T2T-CCS" "HG004-T2T-CCS" "OES103-T2T-10X" "OES103-T2T-CCS" "OES117-T2T-CCS" "OES143-T2T-CCS" "OES148-T2T-CCS" "OES152-T2T-CCS" "OES152-T2T-ONT" "HG003-T2T-CCS" "HT-115-T2T-CCS"  "OES103-T2T-blood" "OES103-T2T-HiC" "OES117-T2T-HiC" "OES143-T2T-HiC" "OES148-T2T-HiC" "OES152-T2T-HiC" )
# output=( "bc1/bc1-ccs" "hg002/hg002-ccs" "hg004/hg004-ccs" "oes103/oes103-10x" "oes103/oes103-ccs" "oes117/oes117-ccs" "oes143/oes143-ccs" "oes148/oes148-ccs" "oes152/oes152-ccs" "oes152/oes152-ont" "hg003/hg003-ccs" "ht115/ht115-ccs" "oes103/oes103-blood" "oes103/oes103-hic" "oes117/oes117-hic" "oes143/oes143-hic" "oes148/oes148-hic" "oes152/oes152-hic" )

files=(  "OES148-T2T-HiC" "OES148-T2T-CCS" )
output=( "oes148-hic" "oes148-ccs" )
for i in "${!files[@]}"; do
	# echo ${files[i]} " goes to "  ${output[i]}
	echo "BEGINNING LAUNCHER ON " ${files[i]}
	./deforest -f ../Coverage/${files[i]}.dat -sigmaResolution 10 -worker 26 -sigmaMin 5 -sigmaMax 25 -o Output/noThin/${output[i]} -L 1000000 -gamma 5 -thin 200 -accelerate 50 -Qmax 6 -alpha 0.1  -smooth 0
	./deforest -f ../Coverage/${files[i]}.dat -sigmaResolution 10 -worker 26 -sigmaMin 5 -sigmaMax 25 -o Output/midThin/${output[i]} -L 1000000 -gamma 5 -thin 200 -accelerate 50 -Qmax 6 -alpha 0.1  -smooth 0.5
	./deforest -f ../Coverage/${files[i]}.dat -sigmaResolution 10 -worker 26 -sigmaMin 5 -sigmaMax 25 -o Output/bigThin/${output[i]} -L 1000000 -gamma 5 -thin 200 -accelerate 50 -Qmax 6 -alpha 0.1  -smooth 0.9
	./deforest -f ../Coverage/${files[i]}.dat -sigmaResolution 10 -worker 26 -sigmaMin 5 -sigmaMax 25 -o Output/sillyThin/${output[i]} -L 1000000 -gamma 5 -thin 200 -accelerate 50 -Qmax 6 -alpha 0.1  -smooth 0.99
	./deforest -f ../Coverage/${files[i]}.dat -sigmaResolution 10 -worker 26 -sigmaMin 5 -sigmaMax 25 -o Output/ludicrousThin/${output[i]} -L 1000000 -gamma 5 -thin 200 -accelerate 50 -Qmax 6 -alpha 0.1  -smooth 0.999
done
