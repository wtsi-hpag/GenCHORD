#!/bin/bash
file=$1

while read p; do
	echo "Launching " $p
  ./deforest -f "DevilData/"$p".dat" -o "Output/"$p"/" -q
done <$file

rm $file