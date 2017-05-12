#!/bin/bash

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified December 5, 2014"
	echo ""
	echo "Description:  Generates random synthetic reads from a reference genome.  Read names indicate their genomic origin."
	echo "Allows precise customization of things like insert size and synthetic mutation type, sizes, and rates."
	echo "Read names generated by this program are used by MakeRocCure (samtoroc.sh) and GradeSamFile (gradesam.sh)."
	echo "They can also be used by BBMap (bbmap.sh) and BBMerge (bbmerge.sh) to automatically calculate"
	echo "true and false positive rates, if the flag 'parsecustom' is used."
	echo ""
	echo "Usage:	randomreads.sh ref=<file> out=<file> length=<number> reads=<number>"
	echo ""
	echo "Basic parameters:"
	echo "out=null        		Output file.  If blank, filename(s) will be autogenerated."
	echo "ref=null        		Reference file.  Not needed if the reference is already indexed."
	echo "build=1         		If multiple references are indexed in the same directory, each needs a unique build ID."
	echo "midpad=300      		Specifies space between scaffolds in packed index."
	echo "reads=0         		Generate this many reads (or pairs)."
	echo "minlength=100        		Generate reads of up to this length."
	echo "maxlength=100        		Generate reads of at least this length."
	echo "length=100        		Generate reads of exactly this length."
	echo "overwrite=true      	 	Set to false to disallow overwriting of existing files."
	echo "replacenoref=false		Set to true to replace N in the reference sequence with random letters."
	echo "illuminanames=false		Set to true to have matching names for paired reads, rather than naming by location."
	echo "prefix=null        		Generated reads will start with this prefix, rather than naming by location."
	echo "seed=0  			Use this to set the random number generator seed; use -1 for a random seed."
	echo ""
	echo "Pairing parameters:"
	echo "paired=false    		Set to true for paired reads."
	echo "interleaved=false		Set to true for interleaved output (rather than in two files)."
	echo "mininsert=      		Controls minimum insert length.  Default depends on read length."
	echo "maxinsert=      		Controls maximum insert length.  Default depends on read length."
	echo "triangle=true   		Make a triangular insert size distribution."
	echo "flat=false      		Make a roughly flat insert size distribution.."
	echo "superflat=false 		Make a perfectly flat insert size distribution."
	echo "gaussian=false  		Make a bell-shaped insert size distribution, with standard deviation of (maxinsert-mininsert)/6."
	echo "samestrand=false		Generate paired reads on the same strand."
	echo ""
	echo "Mutation parameters:"
	echo "snprate=0        		Add snps to reads with this probability (0-1)."
	echo "insrate=0        		Add insertions to reads with this probability (0-1)."
	echo "delrate=0        		Add deletions to reads with this probability (0-1)."
	echo "subrate=0        		Add contiguous substitutions to reads with this probability (0-1)."
	echo "nrate=0          		Add nocalls to reads with this probability (0-1)."
	echo "Note: With a 'rate' of X, each read has an X chance of getting at least 1 mutation, X^2 chance of 2+ mutations,"
	echo "X^3 chance of 3+ mutations, and so forth up to the maximum allowed number of mutations of that type."
	echo ""
	echo "maxsnps=3        		Add at most this many snps per read."
	echo "maxinss=2        		Add at most this many deletions per read."
	echo "maxdels=2        		Add at most this many insertions per read."
	echo "maxsubs=2        		Add at most this many contiguous substitutions per read."
	echo "maxns=0          		Add at most this many blocks of Ns per read."
	echo ""
	echo "maxinslen=12      		Max length of insertions."
	echo "maxdellen=400     		Max length of deletions."
	echo "maxsublen=12      		Max length of contiguous substitutions."
	echo "maxnlen=1        		Min length of N blocks."
	echo ""
	echo "mininslen=1      		Min length of insertions."
	echo "mindellen=1      		Min length of deletions."
	echo "minsublen=2      		Min length of contiguous substitutions."
	echo "minnlen=1        		Min length of N blocks."
	echo ""
	echo "Illumina quality parameters:"
	echo "maxq=36          		Upper bound of quality values."
	echo "midq=32          		Approximate average of quality values."
	echo "minq=28          		Lower bound of quality values."
	echo "q=               		Sets maxq, midq, and minq to the same value."
	echo "adderrors=true    		Add substitution errors based on quality values, after mutations."
	echo ""
	echo "PacBio quality parameters:"
	echo "pacbio=false               	Use a PacBio error model, rather than Illumina error model, and add PacBio errors after mutations."
	echo "pbmin=0.13               	Minimum rate of PacBio errors for a read."
	echo "pbmax=0.17               	Maximum rate of PacBio errors for a read."
	echo ""
	echo "Other Parameters:"
	echo "overlap=1        		Require reads to overlap scaffold end by at least this much."
	echo "randomscaffold=false		Choose random scaffolds without respect to length."
	echo "amp=1            		Simulate highly-amplified MDA single-cell data by setting this to a higher number like 1000."
	echo "replacenoref=false		Replace intra- and inter-scaffold Ns with random bases."
#	echo "colorspace=false  		Generate Solid colorspace reads."
	echo "pbadapter=        		Add adapter sequence to some reads using this literal string."
	echo "fragadapter=      		Add this sequence to paired reads with insert size shorter than read length."
	echo "fragadapter2=     		Use this sequence for read 2."
	echo ""
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       			This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "                		-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx1g"
z2="-Xms1g"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

randomreads() {
	#module unload oracle-jdk
	#module load oracle-jdk/1.7_64bit
	#module load pigz
	local CMD="java $EA $z -cp $CP align2.RandomReads3 build=1 $@"
	echo $CMD >&2
	$CMD
}

randomreads "$@"