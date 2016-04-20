#!/bin/bash

sampleID=$1
fq1=fastq/$sampleID/reads1.fq.gz
fq2=fastq/$sampleID/reads2.fq.gz
#if the fastq files are not trimmed yet.
if [ ! -s fastq/$sampleID/bbduk.stats ]
then
    bbtools=$2/../bbtools # $2 is scriptdir
    bash $bbtools/bbduk.sh in=$fq1 in2=$fq2 ref=$bbtools/full_adaptor_list.fa minlength=36 k=19 ktrim=r hdist=1 mink=12 out=fastq/${sampleID}_trimmed/reads1.fq.gz out2=fastq/${sampleID}_trimmed/reads2.fq.gz stats=fastq/${sampleID}/bbduk.stats overwrite=1
    mv fastq/${sampleID}_trimmed/reads1.fq.gz  $fq1
    mv fastq/${sampleID}_trimmed/reads2.fq.gz $fq2
    rm -rf fastq/${sampleID}_trimmed
fi


