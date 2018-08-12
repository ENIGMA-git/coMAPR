#!/bin/bash
#$ -S /bin/bash

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

config="./config_Sitename.sh"

source $config
## These are all you -- enter your site ID and paths to your files

## how are you running this?? Command-line or on Q-SUB
#Nnodes=${#ROI_LIST[@]} 	# ***uncomment this if using a SGE or PBS cluster *** Set number of nodes to the length of ROI list
Nnodes=1		# *** otherwise we're going to set the number of nodes to 1 and assume you are running locally

#---Set the full path to your R binary
#Rbin=/usr/local/R-3.1.3/bin/R
Rbin="R"
###### optional edits:

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

EXCLUDE_STR="$EXCLUDE_FILE"
#if [ "$EXCLUDE_FILE" != "" ]; then
#	EXCLUDE_STR="-exclude_path $EXCLUDE_FILE"
#else
#	EXCLUDE_STR=""
#fi

#if [ "$METR_PREFIX" != "" ]; then
#	METR_PREFIX_STR="-shape_prefix $METR_PREFIX"
#else
#	METR_PREFIX_STR=""
#fi


#---Section 4. DO NOT EDIT. qsub variable ---
Nroi=${#ROI_LIST[@]}	

SGE_TASK_ID=${PBS_ARRAYID:-$SGE_TASK_ID}

if [ $Nnodes == 1 ]
then
	SGE_TASK_ID=1
fi

if [ ${SGE_TASK_ID} == 1 ]
then
	touch $resDir/roi_list.txt	
	rm $resDir/roi_list.txt
	touch $resDir/roi_list.txt
	roi_text=$(printf "\",\"%s" ${ROI_LIST[@]})
	roi_text=${roi_text:2}"\""

	echo $roi_text >> $resDir/roi_list.txt
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
			${DATADIR}/${METR_PREFIX} \
			${SUBJECTS_COV} \
			${CONFIG_PATH} \
			${dbPath} \
			${METR_LIST} \
			${dbFile} < ${scriptDir}/register_input_csv.R"

	LOG="$logDir/${RUN_ID}_${SITE}_reg_input.log"

	echo $cmd &> $LOG
	echo " " &>> $LOG
	eval $cmd &>> $LOG
	
	
#	LOG="$logDir/${RUN_ID}_${SITE}_run_lm.log"
#	cmd="${Rbin} --no-save --slave --args\
#			${RUN_ID} \
#			${SITE} \
#			${ROI_LIST} \
#			${resDir}  < ${scriptDir}/eregr_run_lm_csv.R"
#
#	echo $cmd &> $LOG
#	echo " " &>> $LOG
#	eval $cmd &>> $LOG
	
done

