#!/bin/bash
#Please change the header here to suit your cluster job submission system
studyName=$1
mkdir -p  $studyName

mkdir -p $studyName/fastq $studyName/resistType $studyName/ErrFiles $studyName/tmpDir $studyName/spadesOutput $studyName/filteredFastq

echo -e $"#!/bin/bash\n#$ -cwd\n#$ -V" >$studyName/run_resistType.sh
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "SCRIPTDIR=$SCRIPTDIR" >>$studyName/run_resistType.sh
cat $SCRIPTDIR/runResistTypeBody.sh >>$studyName/run_resistType.sh
cp $SCRIPTDIR/runBbduk.sh $studyName

echo -e $"#!/bin/bash\n#$ -cwd\n#$ -V" >$studyName/runMLST.sh
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "SCRIPTDIR=$SCRIPTDIR" >>$studyName/runMLST.sh
echo "python $SCRIPTDIR/mlstTyping.py -s $sampleid -c spadesOutput/$sampleid/contigs.fasta" >>$studyName/runMLST.sh

echo -e $"#!/bin/bash\n#$ -cwd\n#$ -V" >$studyName/runISFinder.sh
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "SCRIPTDIR=$SCRIPTDIR" >>$studyName/runISFinder.sh
echo "python $SCRIPTDIR/ISFinder.py -s $sampleid -c spadesOutput/$sampleid/contigs.fasta" >>$studyName/runISFinder.sh

echo -e $"#!/bin/bash\n#$ -cwd\n#$ -V" >$studyName/runPlasmidFinder.sh
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "SCRIPTDIR=$SCRIPTDIR" >>$studyName/runPlasmidFinder.sh
echo "python $SCRIPTDIR/plasmidFinder.py -s $sampleid -c spadesOutput/$sampleid/contigs.fasta" >>$studyName/runPlasmidFinder.sh



cp $SCRIPTDIR/submitJobs.sh $studyName


