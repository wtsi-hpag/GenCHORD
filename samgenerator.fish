#!/bin/fish

cd ../../zn1/stepStone/bam
for f in *.bam
	set base (path change-extension '' $f)
	samtools depth $base.bam | awk '($2%50==0){print $1,$2,$3}' > ../../../jf20/Coverage/$base.dat&
end

