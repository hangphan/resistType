#!/bin/bash
#gradesam in=<infile>

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified May 23, 2014"
	echo ""
	echo "Description:  Creates a ROC curve from a sam file of synthetic reads with headers generated by RandomReads3.java"
	echo ""
	echo "Usage:	samtoroc.sh in=<sam file> reads=<number of reads in input fastq>"
	echo ""
	echo "Parameters:"
	echo "in=<file>     	Specify the input sam file, or stdin."
	echo "thresh=20     	Max deviation from correct location to be considered 'loosely correct'."
	echo "blasr=f       	Set to 't' for BLASR output; fixes extra information added to read names."
	echo "ssaha2=f      	Set to 't' for SSAHA2 or SMALT output; fixes incorrect soft-clipped read locations."
	echo "bitset=t      	Track read ID's to detect secondary alignments."
	echo "              	Necessary for mappers that incorrectly output multiple primary alignments per read."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       	This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "			-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx200m"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
}
calcXmx "$@"

samtoroc() {
	#module unload oracle-jdk
	#module unload samtools
	#module load oracle-jdk/1.7_64bit
	#module load samtools
	local CMD="java $EA $z -cp $CP align2.MakeRocCurve $@"
#	echo $CMD >&2
	$CMD
}

samtoroc "$@"
