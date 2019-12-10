#!/usr/bin/env bash
echo "Usage: Auto-Scaffold.sh <bam file> <alignment> <optional new prefix>"
in_bam=$1
ref=$2
prefix=$3
t1=$(mktemp)
t2=$(mktemp)
t3=$(mktemp)
trap "rm -f ${t1} ${t2} ${t3} " 0 2 3 15 ## remove temp files when we die
if [[ $prefix == "" ]]   # if the user wants a new prefix put it in if not use the name from the bam file
 then prefix=${in_bam::-3}
fi

echo "Calculating coverage"
bamCoverage -b ${in_bam} -of bedgraph -bs 1 -o ${t1} & ## make a coverage map as a bedgraph
echo "Calling variants"
freebayes -C 1 -f ${ref} ${in_bam} > ${prefix}.vcf  ## call variants form ref
bgzip ${prefix}.vcf
bcftools index ${prefix}.vcf.gz
echo "Making consensus and masking"
awk '{if ($4==0){print $0;}}' ${t1} > ${t2}  # convert to only region we need to mask later
bcftools consensus -f ${ref} -s unknown  ${prefix}.vcf.gz > ${t3} ## make an unmasked consensus from the vcf
bedtools maskfasta -fi ${t3} -bed ${t2} -fo ${prefix}_Scaffold.fa ## mask regions that had zero coverage
