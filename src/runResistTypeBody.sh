

while [[ $# > 1 ]]
do
key="$1"

case $key in
    -d|--bamdir)
    BAMDIR="$2"
    shift # past argument
    ;;
    -s|--sampleid)
    SAMPLEID="$2"
    shift # past argument
    ;;
    -t|--nthread)
    NTHREAD="$2"
    shift # past argument
    ;;

    -r|--refid)
    REFID="$2" #if your sample is Ecol or Kpne, please indicate this using -r flag
    shift # past argument
    ;;

    -p|--plateid)
    PLATEID="$2"
    shift # past argument
    ;;
    -b|--bamFile)
    BAMFILE="$2"
    shift # past argument
    ;;
    -m|--metagenomics)
    METAG="$2"
    shift # past argument
    ;;
    -c|--contig)
    CONTIG="$2"
    shift # past argument
    ;;

    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

if [ $PLATEID ]
then
    SAMPLEID=$PLATEID/$SAMPLEID
fi

if [ ! $NTHREAD ]
then
    NTHREAD=1
fi

if [ $BAMDIR  ]
then    
    #if provided with a directory containing bam files to process
    #assuming of this form data/OX1173_ECOL_APHA_Molsig/raw_input/MG1655_ref
    #this only work for data from MMM pipeline system
    #IFS="/" read -a array <<< "$BAMDIR"
    #echo ${array[3]}
    #SAMPLEID=${array[${#array[@]}-1]}
    #PLATEID=${array[1]}
    #SAMPLEID=$PLATEID/$SAMPLEID
    #echo $SAMPLEID
    shopt -s nullglob
    array=($BAMDIR/*.bam)
    bamFiles=$(printf " %s" "${array[@]}")
    echo ${#array[@]} 
    if [ ${#array[@]} -gt 1 ]
    then
	samtools merge tmpDir/$SAMPLEID.bam $bamFiles
	BAMFILE=tmpDir/$SAMPLEID.bam
    else
	BAMFILE=${array[0]}
    fi
    echo $BAMFILE

fi

if [ $BAMFILE ] 
then
#    only use this with the MMM pipeline system
#    IFS="/" read -a array <<< "$BAMFILE"
#    bamFileName=${array[${#array[@]}-1]}
#    SAMPLEID="${bamFileName%_R*}"
#    REFID=${array[${#array[@]}-3]}
    echo Do nothing
fi

if [ ! $SAMPLEID ]
then
    echo "No sample id detected. Quit here"
    exit
fi

for f in fastq spadesOutput tmpDir resistType filteredFastq
do
    mkdir -p $f/$SAMPLEID
done

fq1=fastq/${SAMPLEID}/reads1.fq.gz
fq2=fastq/${SAMPLEID}/reads2.fq.gz

if [ ! -s $fq2 ]
then
    echo No fastq file detected. Try to see if bam file is provided to generate fastq from here.
    if [ ! -s $BAMFILE ]
    then
	echo No bam file provided to generate reads. Trying to look for assembly. 
    else
	if [ ! -s $CONTIG ]
	then
	    echo No assembly file provided. Quit here 
	    exit
	fi
    fi

    if [ $BAMFILE ]
    then
	echo Running samtools 
	$S

CRIPTDIR/../bin/samtools sort -n $BAMFILE tmpDir/$SAMPLEID.temp 
	$SCRIPTDIR/../bin/bedtools bamtofastq -i tmpDir/$SAMPLEID.temp.bam -fq fastq/$SAMPLEID/reads1.fq -fq2 fastq/$SAMPLEID/reads2.fq  
	gzip -f  fastq/$SAMPLEID/reads1.fq 
	gzip -f  fastq/$SAMPLEID/reads2.fq 
	rm  tmpDir/$SAMPLEID.temp.bam  
    fi

fi

#trimming read adaptors if not trimmed yet
sh runBbduk.sh $SAMPLEID $SCRIPTDIR

if [ ! $REFID ]
then 
    REFID=Ente #fake refid
fi

rm resistType/$SAMPLEID/ref*

if [ ! $METAG ] #not running metagenomic based typing
then
    if [ ! -s spadesOutput/$SAMPLEID/contigs.fasta ]
    then
	spadesPath=$SCRIPTDIR/../bin/spades.py
	spadesPath=spades.py
	outDir=tmpDir/$SAMPLEID
	$spadesPath -1 fastq/$SAMPLEID/reads1.fq.gz -2 fastq/$SAMPLEID/reads2.fq.gz -o $outDir   --careful -t $NTHREAD
	mv $outDir/*.fas* spadesOutput/$SAMPLEID
	mv $outDir/scaffold* spadesOutput/$SAMPLEID
	mv $outDir/params.txt spadesOutput/$SAMPLEID
	mv $outDir/spades.log spadesOutput/$SAMPLEID
	rm -rf $outDir
    fi
    if [ -s spadesOutput/$SAMPLEID/contigs.fasta ]
    then
	python $SCRIPTDIR/resistType_v0.1.py -s $SAMPLEID -r $REFID
	sh runISFinder.sh $SAMPLEID
	sh runPlasmidFinder.sh $SAMPLEID
	sh runMLST.sh $SAMPLEID
    else
	echo SPAdes assembly failed, running resistance prediction on metagenomics mode no copy number estimation
	python $SCRIPTDIR/resistType_v0.1.py -s $SAMPLEID -m 1
	sh runMLST.sh $SAMPLEID
	sh runISFinder.sh $SAMPLEID
	sh runPlasmidFinder.sh $SAMPLEID
    fi
else
    python $SCRIPTDIR/resistType_v0.1.py -s $SAMPLEID -m 1
fi


