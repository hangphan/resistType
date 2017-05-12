#!/bin/bash
#grademerge in=<infile> out=<outfile>

function usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2015

Description:  Grades correctness of merging synthetic reads with headers 
              generated by RandomReads and re-headered by RenameReads.

Usage:  grademerge.sh in=<file>

Parameters:
in=<file>           Specify the input file, or 'stdin'.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

function grademerge() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA -Xmx200m -cp $CP jgi.GradeMergedReads $@"
#	echo $CMD >&2
	$CMD
}

grademerge "$@"
