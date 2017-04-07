#Registering input data script for mass_uv_regr_v20 package
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

suppressPackageStartupMessages(source(file='lib_eregr_core.R'))

cmdargs = commandArgs(trailingOnly=T)
if(length(cmdargs) != 11) {
	dput(cmdargs)
	stop("incorrect number of parameters")
}

study_Id <- cmdargs[1]

site_Id <- cmdargs[2]

metr_folder <- cmdargs[3]

logDir <- cmdargs[4]

resDir <- cmdargs[5]

subjects_cov <- cmdargs[6]

ROI_arg <- cmdargs[7]

config_path <- cmdargs[8]   #docs.google 

Exclude_Path <- cmdargs[9]

QA_LEVEL <- as.numeric(cmdargs[10])

db_path <- cmdargs[11]


LOG_FILE<-paste(logDir, '/',RUN_ID,'_',SITE,'.log',sep='')
Results_Path <- paste(resDir,'/',RUN_ID,'_',sep='')

message("1. Session Info")
sessionInfo()


# Connect to db
message("2. Connecting to database")
#now it's configured for only one database, I'll fix that later
con <- eregr_connect("db_mysql",RMySQL::MySQL(),user='eregr',password='eregr',dbname='apr7demo',host='127.0.0.1',port=3306)
#con <- eregr_connect(db_path)

message("3. Registering study")
is_study_reg <- tryCatch(eregr_register_study(con,study_Id = study_Id, study_Name = '',gsheet_path = config_path),
	error = function(e) {
		message("Message: Cannot register study ", study_Id, " ", e$message,"\n")
		FALSE
		}
	)
if(is_study_reg)
	message("Study ", study_Id, " registered succesfully")
# here we can check the existence of study

message("4. Registering site")
is_site_reg <- tryCatch(eregr_register_site(con, study_Id = study_Id, site_Id = site_Id),
	error = function(e) {
		message("Message: Cannot register site ", site_Id, " ", e$message)
		FALSE
		}
)
if(is_site_reg)
	message("Site ", site_Id, " registered successfully")

message("5. Registering covariates")
cov_data <- tryCatch(eregr_register_covariates(con,study_Id,site_Id,subjects_cov),
		     error = function (e) {
				message("Error: could not register covariates")
				stop(e)
		     })
message("Covariates registered successfully")

message ("6. Registering metrics files")

ROI_vec <- str_split(ROI_arg,",")[[1]]
subj_list <- eregr_get_subjlist_from_covars(con,study_Id,site_Id)
study_metr <- eregr_get_study_metrics(con,study_Id)
is_metr_read <- tryCatch(eregr_read_shape_all_subjects(con,metr_folder, subj_list$subjID,
                        as.character(ROI_vec),study_metr$metr_name,study_Id,site_Id,Exclude_Path),
		   error = function(e) {
			message("Error: could not register metric files")
			stop(e)
	
		})
message("Metrics files registered successfully")
message("7. Disconnecting from database")

is_disconnected <- eregr_disconnect(con)

message("Disconnected succesfully")
message("Finished successfully")
