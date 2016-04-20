#!/bin/bash
#Please change the header here to suit your cluster job submission system
studyName=$1
mkdir -p  $studyName

mkdir -p $studyName/fastq $studyName/resistType $studyName/ErrFiles $studyName/tmpDir $studyName/spadesOutput $studyName/filteredFastq
echo "#!/bin/bash" >$studyName/run_resistType.sh
echo "#$ -cwd" >>$studyName/run_resistType.sh
echo "#$ -V" >>$studyName/run_resistType.sh
echo "#$ -P bag.prjb -q short.qb ">>$studyName/run_resistType.sh
echo "#$ -pe shmem 1 " >>$studyName/run_resistType.sh
echo "#$ -e $PWD/ErrFiles/" >>$studyName/run_resistType.sh
echo "#$ -o $PWD/ErrFiles/" >>$studyName/run_resistType.sh
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "SCRIPTDIR=$SCRIPTDIR" >>$studyName/run_resistType.sh
cat $SCRIPTDIR/runResistTypeBody.sh >>$studyName/run_resistType.sh
cp $SCRIPTDIR/runBbduk.sh $studyName
cp $SCRIPTDIR/submitJobs.sh $studyName


