#---DO NOT CHANGE. SCRIPT PATH
scriptDir="/path_to_scripts" ## where you have downloaded the ENIGMA Regression R scripts!
dbPath="NA"
dbFile="database_connect_YOURSTUDY.R"



#---Study level config
## Get the following from your working group leader ## 
RUN_ID="STUDY_ID"
CONFIG_PATH="https://docs.google.com/spreadsheets/d/1_5937FsjsSp7AEYkXQ3i67VBJx6Kq9b3qHyq7xXVPm0"
ROI_LIST="R1,R2,R3,.." # e.g. for DTI ANALYSIS: "ACR,ACR_L,ACR_R,ALIC,ALIC_L,ALIC_R,AverageFA,BCC,CC,CGC,CGC_L,CGC_R,CGH,CGH_L,CGH_R,CR,CR_L,CR_R,CST,CST_L,CST_R,EC,EC_L,EC_R,FX,FX_ST_L,FX_ST_R,FXST,GCC,IC,IC_L,IC_R,IFO,IFO_L,IFO_R,PCR,PCR_L,PCR_R,PLIC,PLIC_L,PLIC_R,PTR,PTR_L,PTR_R,RLIC,RLIC_L,RLIC_R,SCC,SCR,SCR_L,SCR_R,SFO,SFO_L,SFO_R,SLF,SLF_L,SLF_R,SS,SS_L,SS_R,UNC,UNC_L,UNC_R"


#Site level config
logDir="/path_to_site/logs"      ## directory to be created to output the log files
resDir="/path_to_site/res"      ## directory to be created to output the log files
SITE="Sitename"
DATADIR="/path_to_site/data"
SUBJECTS_COV="/path_to_site/data/cov.csv"
METR_LIST="M1,M2,M3,M4" # e.g "FA,MD,RD,AD"
METR_PREFIX="metr_"


