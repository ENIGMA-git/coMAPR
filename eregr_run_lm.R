#Running linear models and writing results to database script for mass_uv_regr_v20 package
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
if(length(cmdargs) != 7) {
	dput(cmdargs)
	stop("incorrect number of parameters")
}

study_Id <- cmdargs[1]

site_Id <- cmdargs[2]

CURRENT_ROI <- cmdargs[3]

Exclude_Path <- cmdargs[4]

QA_LEVEL <- cmdargs[5]

Results_Path <- cmdargs[6]

db_path <- cmdargs[7]

message("1. Session Info")
sessionInfo()


# Connect to db
message("2. Connecting to database")
#con <- eregr_connect(db_path)

# now it's configured for only one database, I'll fix that later
con <- eregr_connect("db_mysql",RMySQL::MySQL(),user='eregr',password='eregr',dbname='apr7demo',host='127.0.0.1',port=3306)



# check that study_exist
# check that site exist


message("3. Reading QC file")
exclude_csv <- tryCatch(read.csv(Exclude_Path,stringsAsFactors = FALSE),
			error = function(e) {
				message("Error: Cannot read QC file")
				stop(e)
			}
			)
str_search <- paste("SubjID|ROI",CURRENT_ROI,sep="")
str_filter <- paste("ROI",CURRENT_ROI,"<",as.character(QA_LEVEL),sep="")
exclude_idx <- which(str_detect(names(exclude_csv),str_search))
if(length(exclude_idx) !=2) stop ("could not filter ROI and SubjID from QC file.")
# here we want to be realy carefull and exit even on warning:
options(warn=2)
subjectsExcl <- exclude_csv %>% 
			select(exclude_idx) %>%
				filter_(str_filter) %>%
					rename(subjID = SubjID) %>%
						select(subjID)
#switching warning priority back:
options(warn=1)					
message ("QC file read succesfully")

message ("4. Running linear models")
tic()
df_list <- NULL 
lm_err <- try(df_list <- eregr_run_linear_models(con, study_Id, site_Id, CURRENT_ROI,subjectsExcl))
toc()
if(is.null(df_list))
	stop(lm_err)

message("4.1 Model computation finished. Brief summary:")
errs <- map(df_list,"error")
for (i in seq_along(errs)) {
    if (is.null(errs[[i]])) 
        message("Success: \t", names(errs)[[i]])
    else 
        message("Failed: \t", names(errs)[[i]])
}

message("5. Writing model results to database")

res <- tryCatch(eregr_write_results_linear_models(con,study_Id,site_Id,CURRENT_ROI,df_list),
	error = function (e) {
			message("Error: couldn't write results to database")
			stop(e)

		})
message("6. Disconnecting from database")
is_disconnected <- eregr_disconnect(con)

message("7. Saving errors, models and results to files")
saveRDS(errs,paste("error_",study_Id,"_",CURRENT_ROI,"_",site_Id,".rds",sep=""))

# if needed we can save all models for debug purposes
#saveRDS(df_list,paste("lm_",study_Id,"_",CURRENT_ROI,"_",site_Id,".rds",sep=""))

saveRDS(res,paste("res_",study_Id,"_",CURRENT_ROI,"_",site_Id,".rds",sep=""))

message("Finished successfully")
