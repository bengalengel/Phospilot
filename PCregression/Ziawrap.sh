#!/usr/bin/bash

####below are commands for the SGE

#this is an array job with "tasks" numbered 1 through 5
#### -t 1-5

# the qstat job name
#$ -N RTest

# use the real bash shell
#$ -S /bin/bash

# use my environment variables
#$ -V

# Send the output to the directory the bashwrapper.sh file is in
#$ -cwd

#concatenate the std out and in 
#$ -j y

# mail me ...Doesn't work!
#$ -M bengelmann@uchicago.edu

# ... when the job (b)egins, (e)nds (a)borts, or (s)uspends
#$ -m beas

#####Call the R script without invoking 'module load' because I have passed qsub "V".

Rscript $HOME/for_brett/ZiaDataPCregress.R 12




