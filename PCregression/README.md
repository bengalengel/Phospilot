These scripts perform PCregression/pQTL optimization on the GelPrep protein group estimates.
The scripts are called from "NormProt.R" when the required files/folders are missing.
"PC_regress" is to be run with Rscript using the "ZiaArrayWrap.sh" wrapper script on the cluster as a batch array job.
The results are manually moved to the "$HOME/Phospilot/PCregress/Output" directory.
If not working with cluster mounted, this output should be on the drive performing the analysis.