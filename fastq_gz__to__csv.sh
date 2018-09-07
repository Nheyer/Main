#!/usr/bin/env bash

##### This file should be placed in the directory with all the gziped R1 fastq files you want with R2 in a seperate directory named the same
##### The directures should both b in the same directory
##### the path to the reference you want to use then threads
##### (should be renamed ref.fasta) and the python script 'temps_to_csv_line.py'
##### The user should append the output of this file to a file

# define the 3 different methods people were using
no_buf=' chr1:119911553-120069703 chr1:146151908-146229032 chr1:148602850-148679774 chr1:149390591-149468389 chr1:120723917-120798050'
mid_buf=' chr1:119989248-120190000 chr1:149328818-149471561 chr1:148600079-148801427 chr1:146149145-146328264 chr1:120705669-120801220'
big_buf='chr1:143200001-150600000 chr1:120400001-121700000'
#######~~~~~~~~~~~~ FOR THE USER TO EDIT ~~~~~~~~~~~~~#######
regions=$no_buf
picard="/public/home/nheyer/bin/picard/build/libs/picard.jar"
ref=$1
t=$2
#### startup ####
echo 'Sample , reads , duplicates , 1q21.1 , 1q21.1 Enrichment , % mapped , %duplicates , unique BCs , reads/unique BC , average insert size , path to histogram '
# edit this line to put the absolute path to the directory then *.fastq.gz
for file in *.fastq.gz
 # unzip the files
 do pigz -p $t -d $file
 R2=${file%%_R1_*}_R2_${file##*_R1_} 
 pigz -p $t -d ../R2/$R2
 #### processing ####
 bwa mem -t $t $ref "${file::-3}"  "../R2/${R2::-3}" > temp_1-2.sam
 samtools view -Sb temp_1-2.sam > temp_2-2.bam
 samtools sort temp_2-2.bam -o temp_3-2.bam
 # use the original file name (without the endings) to name the file and give it a .infoa file ending
 samtools flagstat temp_3-2.bam > "${file::-9}.infoa"
 # create an index then take only the reads that are in the defined regions and put it into a new temp file, then resort and put it into yet another temp file
 samtools index temp_3-2.bam
 samtools view -b temp_3-2.bam $regions > temp_4.bam
 samtools sort temp_4.bam -o temp_5.bam
 java -jar $picard CollectInsertSizeMetrics I=temp_5.bam O="${file::-9}.info_insert" H="${file::-9}_histogram.pdf"
 # use the original file name (without the endings) to name the flagstat output from data in the bam containing only the regions defined above and give it a .infob file ending
 samtools flagstat temp_5.bam > "${file::-9}.infob"
 # run the python script to take the info we need and convert it to one coma seperated line
 python3 temps_to_csv_line_w_sets.py -ia "${file::-9}.infoa" -ib "${file::-9}.infob" -iq "${file::-3}" -ii "${file::-9}.info_insert"


 #### cleaning ####
 # tell the user what we are doing using the standard error
 echo "cleaning up" 1>&2
 # remove all the files we made and re zip the fastq file
 rm temp_*
 pigz -p $t "${file::-3}"
 pigz -p $t "../R2/${file::-3}"
 done
