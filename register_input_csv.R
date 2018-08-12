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
if(length(cmdargs) != 8) {
	dput(cmdargs)
	stop("incorrect number of parameters")
}

study_Id <- cmdargs[1]

site_Id <- cmdargs[2]

metr_folder <- cmdargs[3]

subjects_cov <- cmdargs[4]

config_path <- cmdargs[5]   #docs.google 

db_path <- cmdargs[6]

metr_arg <- cmdargs[7]

db_connect_path <- cmdargs[8]

message("1. Session Info")
sessionInfo()

# Connect to db
message("2. Connecting to database")
source(file=db_connect_path)
con <- database_connect()

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
metr_list <- str_split(metr_arg,",")[[1]]
for (i in seq_along(metr_list)) {
	cur_metr_fname <- paste(metr_folder,metr_list[[i]],'.csv',sep='')

	is_metr_read <- tryCatch(eregr_register_metrics_csv(con,study_Id,site_Id,cur_metr_fname,metr_list[[i]]),
	  		 error = function(e) {
				message("Error: could not register metric files")
				stop(e)
			})

}


message("Metrics files registered successfully")

# here add feature set registration and new regressors computation

message("7. Registering feature sets")
fs_data<- tryCatch (eregr_register_feature_sets(con,study_Id,site_Id),
			error = function(e) {
				message("Error: could not register feature sets")
				stop(e)
			})
message("Feature sets registered successfully")

message("8. Computing new covariates for feature sets")
new_covar <- tryCatch(eregr_int_compute_new_covariates(con,fs_data,study_Id,site_Id),
			error=function(e) {
				message("Error: could not compute new covariates for feature sets")
				stop(e)
			})


message("9. Disconnecting from database")

is_disconnected <- eregr_disconnect(con)

message("Disconnected succesfully")
message("Finished successfully")
