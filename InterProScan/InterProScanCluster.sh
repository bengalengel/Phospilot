#!/usr/bin/bash

##This wrapper script it submitted from the head node WITHIN interproscan directory. qsub -l h_vmem=6g InterProScanCluster.sh

####below are commands for the SGE

#this is an array job with "tasks" numbered 1 through N (for the N fasta files)
#$ -t 1-5

# the qstat job name
#$ -N IPSArraytest

# use the real bash shell
#$ -S /bin/bash

# use my environment variables
#$ -V

# Set the working directory to the 'for_brett' working directory
#$ -wd /mnt/lustre/home/bengelmann/Programs/my_interproscan/interproscan-5.14-53.0

#concatenate the std out and in 
#$ -j y

# mail me ...Doesn't work!
#$ -M bengelmann@uchicago.edu

# mail me when the job (b)egins, (e)nds (a)borts, or (s)uspends
#$ -m beas


#####Call interproscan.sh on each fasta file ot identify pfam domains and write .tsv files to same directory if they are not present

./interproscan.sh -i ./FASTA/Interpro.$SGE_TASK_ID.fa -appl Pfam -f tsv -o ./FASTA/Interpro.out.$SGE_TASK_ID.tsv


