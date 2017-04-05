##!/bin/bash
#$ -S /bin/bash
#$ -o /ifshome/disaev/ENIGMA_TUTORIAL/log/main_log.log 

#----Wrapper for csv version of mass_uv_regr.R script.
#----See readme to change the parameters for your data
#-Dmitry Isaev
#-Boris Gutman
#-Neda Jahanshad
# Beta version for testing on sites.
#-Imaging Genetics Center, Keck School of Medicine, University of Southern California
#-ENIGMA Project, 2015
# enigma@ini.usc.edu 
# http://enigma.ini.usc.edu
#-----------------------------------------------

#---Section 1. Script directories
scriptDir="/ifshome/disaev/projects/mass_uv_regr_v20/gitrepo" ## where you have downloaded the ENIGMA Regression R scripts!
resDir="/ifs/loni/faculty/thompson/four_d/disaev/projects/mass_uv_regr_20/enigma_schizo_shape_v20/res"      ## directory to be created for your results!
logDir="/ifs/loni/faculty/thompson/four_d/disaev/projects/mass_uv_regr_20/enigma_schizo_shape_v20/logs"        ## directory to be created to output the log files
dbPath="/ifs/loni/faculty/thompson/four_d/disaev/projects/mass_uv_regr_20/sandbox/teststand.sqlite"

#---Section 2. Configuration variables-----
## Get the following from your working group leader ## 
RUN_ID="TESTSTAND_v2"
ROI_LIST=("10" "58")
############
## These are all you -- enter your site ID and paths to your files
SITE="SantMonica"

## how are you running this?? Command-line or on Q-SUB
#Nnodes=${#ROI_LIST[@]} 	# ***uncomment this if using a SGE or PBS cluster *** Set number of nodes to the length of ROI list
Nnodes=1		# *** otherwise we're going to set the number of nodes to 1 and assume you are running locally

#---Set the full path to your R binary
#Rbin=/usr/local/R-3.1.3/bin/R
Rbin="R"
###### optional edits:
QA_LEVEL="2"
EXCLUDE_FILE="/ifs/loni/faculty/thompson/four_d/Artemis/Shape/ENIGMA_SZ/Galway/Galway_output/QA_status_unix.csv"

######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## no need to edit below this line ##########
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
#---Section 3. DO NOT EDIT. some additional processing of arbitrary variables

EXCLUDE_STR="$EXCLUDE_FILE"

#---Section 4. DO NOT EDIT. qsub variable ---
Nroi=${#ROI_LIST[@]}	

SGE_TASK_ID=${PBS_ARRAYID:-$SGE_TASK_ID}

if [ $Nnodes == 1 ]
then
	SGE_TASK_ID=1
fi

if [ ${SGE_TASK_ID} == 1 ]
then
	touch $scriptDir/roi_list.txt	
	rm $scriptDir/roi_list.txt
	touch $scriptDir/roi_list.txt
	roi_text=$(printf "\",\"%s" ${ROI_LIST[@]})
	roi_text=${roi_text:2}"\""

	echo $roi_text >> $scriptDir/roi_list.txt
fi
NchunksPerTask=$((Nroi/Nnodes))
start_pt=$(($((${SGE_TASK_ID}-1))*${NchunksPerTask}+1))
end_pt=$((${SGE_TASK_ID}*${NchunksPerTask}))

if [ "$SGE_TASK_ID" == "$Nnodes" ]
then
end_pt=$((${Nroi}))
fi

#---Section 6. DO NOT EDIT. Running the R script
cd ${scriptDir}
for ((i=${start_pt}; i<=${end_pt};i++));
do
	cur_roi=${ROI_LIST[$i-1]}  
	cmd="${Rbin} --no-save --slave --args\
			${RUN_ID}\
			${SITE} \
			${cur_roi} \
			${EXCLUDE_STR} \
			${QA_LEVEL} \ 
			${resDir} \
			${dbPath} < ${scriptDir}/eregr_run_lm.R"

	LOG="$logDir/${RUN_ID}_${cur_roi}_${SITE}_eregr_lm.log"

	echo $cmd &> $LOG
	echo " " &>> $LOG
	eval $cmd &>> $LOG
done

