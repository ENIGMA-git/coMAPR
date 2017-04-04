#Main script for mass_uv_regr_v20 package
# --Dmitry Isaev--
# --Boris Gutman--
# --Neda Jahanshad--
# --Imaging Genetics Center, Keck School of Medicine, University of Soutrhern California
# --ENIGMA Project, 2017
# enigma@ini.usc.edu
# http://enigma.ini.usc.edu
#
# beta version being tested now
# report bugs to GitHub issue tracker for this project

source(file='lib_eregr_core.R')

cmdargs = commandArgs(trailingOnly=T)
ID <- cmdargs[1]
#RUN_ID=cmdargs[1]

SITE <- cmdargs[2]

DATADIR <- cmdargs[3]

logDir <- cmdargs[4]
LOG_FILE<-paste(logDir, '/',RUN_ID,'_',SITE,'.log',sep='')

resDir <- cmdargs[5]
Results_CSV_Path <- paste(resDir,'/',RUN_ID,'_',sep='')

subjects_cov <- cmdargs[6]
Subjects_Path <- subjects_cov

CURRENT_ROI <- cmdargs[7]
RUN_ID <- paste(ID,'_',CURRENT_ROI,sep='')

Config_Path <- cmdargs[8]   #docs.google 

dp_path <- cmdargs[9]

#db_path="/ifs/loni/faculty/thompson/four_d/disaev/projects/mass_uv_regr_20/sandbox/teststand.sqlite"

con <- tryCatch(eregr_connect(db_path),
		error = function (e) stop(e))

eregr_disconnect(con)
