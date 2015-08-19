#!/usr/bin/bash

####below are commands for the SGE

#this is an array job with "tasks" numbered 1 through 21 (for 20 pcs)
#$ -t 1-31

# the qstat job name
#$ -N PCregress

# use the real bash shell
#$ -S /bin/bash

# use my environment variables
#$ -V

# Set the working directory to the 'for_brett' working directory
#$ -wd /mnt/lustre/home/bengelmann/Phospilot/PCregression/

#concatenate the std out and in 
#$ -j y

# mail me ...Doesn't work!
#$ -M bengelmann@uchicago.edu

# mail me when the job (b)egins, (e)nds (a)borts, or (s)uspends
#$ -m beas

##edit so I can pass task '1' as '0' PCs regressed
i=$(expr $SGE_TASK_ID - 1)

#####Call the R script without invoking 'module load' because I have passed qsub "V". Pass the number of PCs to regress away from protein quants

Rscript $HOME/Phospilot/PCregression/PC_regress.R $i




