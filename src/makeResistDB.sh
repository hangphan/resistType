
#preprocessing the database, generating fasta files, forward DNA sequences, protein sequences, convert gene names to computer friendly ones
python src/preprocessDB.py $1
#cluster genes using cd hit
bin/cd-hit-est -i resistDB/ResistanceGeneSeqs_fw.fasta -o resistDB/temp90.fa -c 0.90 -n 8 -aS 0.8 -g 1 -p 1 -d 0 >resistDB/db_90.log
#generate fasta files for each gene clusters found in the clustering step
rm resistDB/clusterFas/*
mkdir -p resistDB/clusterFas/
python src/makeClusterMFA.py
#prepare resistance gene file for bwa mapping
bwa index resistDB/ResistanceGeneSeqs.fasta
exit
