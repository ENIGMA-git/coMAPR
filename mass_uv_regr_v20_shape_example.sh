#!/bin/bash
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
scriptDir="/ifshome/disaev/ENIGMA_TUTORIAL/scripts/" ## where you have downloaded the ENIGMA Regression R scripts!
resDir="/ifshome/disaev/ENIGMA_TUTORIAL/results/"   ## directory to be created for your results!
logDir="/ifshome/disaev/ENIGMA_TUTORIAL/log/"        ## directory to be created to output the log files

#---Section 2. Configuration variables-----
## Get the following from your working group leader ## 
RUN_ID="TESTSTAND_v2"
CONFIG_PATH="https://docs.google.com/spreadsheets/d/1-ThyEvz1qMOlEOrm2yM86rD_KABr_YE4yqYmHogaQg0/pubhtml"
ROI_LIST=("10" "58")
############

## These are all you -- enter your site ID and paths to your files
SITE="MDR"
DATADIR="/ifs/loni/faculty/thompson/four_d/Artemis/Shape/ENIGMA_SZ/Galway/Galway_output"
SUBJECTS_COV="/ifs/loni/faculty/thompson/four_d/Artemis/Shape/ENIGMA_SZ/Galway/Galway_output/Covariates_vols.csv"

## how are you running this?? Command-line or on Q-SUB
Nnodes=${#ROI_LIST[@]} 	# ***uncomment this if using a SGE or PBS cluster *** Set number of nodes to the length of ROI list
#Nnodes=1		# *** otherwise we're going to set the number of nodes to 1 and assume you are running locally

#---Set the full path to your R binary
Rbin=/usr/local/R-3.1.3/bin/R

###### optional edits:
QA_LEVEL=""
EXCLUDE_FILE="/ifs/loni/faculty/thompson/four_d/Artemis/Shape/ENIGMA_SZ/Galway/Galway_output/QA_status_unix.csv"
METR_PREFIX="metr_"

######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## no need to edit below this line ##########
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
#---Section 3. DO NOT EDIT. some additional processing of arbitrary variables

if [ ! -d $scriptDir ]
then
   mkdir -p $scriptDir
fi

if [ ! -d $resDir ]
then
   mkdir -p $resDir
fi

if [ ! -d $logDir ]
then
   mkdir -p $logDir
fi

if [ "$EXCLUDE_FILE" != "" ]; then
	EXCLUDE_STR="-exclude_path $EXCLUDE_FILE"
else
	EXCLUDE_STR=""
fi

if [ "$METR_PREFIX" != "" ]; then
	METR_PREFIX_STR="-shape_prefix $METR_PREFIX"
else
	METR_PREFIX_STR=""
fi


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

OUT=$scriptDir/log.txt
touch $OUT
for ((i=${start_pt}; i<=${end_pt};i++));
do
	cur_roi=${ROI_LIST[$i-1]}  
	cmd="${Rbin} --no-save --slave --args\
			${RUN_ID}\
			${SITE} \
			${DATADIR} \
			${logDir} \
			${resDir}
			${SUBJECTS_COV} \
			${cur_roi} \
			${CONFIG_PATH} \
			${QA_LEVEL}  \
			${EXCLUDE_STR} \
			${METR_PREFIX_STR} \
			<  ${scriptDir}/mass_uv_regr.R
		"
	echo $cmd
	echo $cmd >> $OUT
	eval $cmd
done
