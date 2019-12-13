#!/usr/bin/env bash
echo "Usage: Auto-Scaffold.sh <fasta> <alignment> <optional new prefix>"
in_fasta=$1
ref=$2
prefix=$3
t1=$(mktemp)
t2=$(mktemp)
t3=$(mktemp)
trap "rm -f ${t1} ${t2} ${t3} " 0 2 3 15 ## remove temp files when we die
if [[ $prefix == "" ]]   # if the user wants a new prefix put it in if not use the name from the bam file
 then prefix=${in_fasta::-2}
fi
in_bam=${prefix}sorted.bam

if [$(grep -v \> $ref | fold -w1 | grep [^ATCGNatcgn] -c) > 0]; then
 echo "We have a non-ATCGN base in the reference, we will make a temporary ref to deal with this, if there are more then 1 strand in ref it will fail"
 n_ref=$(mktemp)
 trap "rm -f ${n_ref}"
 grep \> ${ref} > ${n_ref}
 cat ${ref} | grep -v ">" | sed 's/[^ATCGNatcgn]/N/g' >> ${n_ref}
 ref=$n_ref
fi
echo"Aligning to the reference"
bwa mem $ref $in_fasta > ${prefix}sam
samtools view -Sb ${prefix}sam > ${prefix}bam
samtools sort ${prefix}bam -o $in_bam
samtools index $in_bam
echo "Calculating coverage"
bamCoverage -b ${in_bam} -of bedgraph -bs 1 -o ${t1} & ## make a coverage map as a bedgraph
echo "Calling variants"
freebayes -C 1 -f ${ref} ${in_bam} > ${prefix}.vcf  ## call variants form ref
bgzip ${prefix}.vcf
bcftools index -f ${prefix}.vcf.gz
echo "Making consensus and masking"
awk '{if ($4==0){print $0;}}' ${t1} > ${t2}  # convert to only region we need to mask later
bcftools consensus -f ${ref} -s unknown  ${prefix}.vcf.gz -o ${t3} ## make an unmasked consensus from the vcf
bedtools maskfasta -fi ${t3} -bed ${t2} -fo ${prefix}_Scaffold.fa ## mask regions that had zero coverage
