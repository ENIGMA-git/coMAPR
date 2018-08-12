#Running meta-analysis and writing results to database script for mass_uv_regr_v20 package
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

site_VEC_arg <- cmdargs[2]
site_VEC <- str_split(site_VEC_arg,",")[[1]]

ROI_VEC_arg <- cmdargs[3]
ROI_VEC <- str_split(ROI_VEC_arg,",")[[1]]
names(ROI_VEC) <- ROI_VEC

#CURRENT_ROI <- cmdargs[3]

Results_Path <- cmdargs[4]

dbFile <- cmdargs[5]

n_param <- cmdargs[6]

n_min <- as.numeric(cmdargs[7])

subjectsExcl=NA

message("1. Session Info")
sessionInfo()


# Connect to db
message("2. Connecting to database")
source(file=dbFile)
con <- database_connect()

print(study_Id)
print(site_VEC)
print(ROI_VEC)
print(Results_Path)

message("3. Getting active linear models to meta-analyze")

study_info <- eregr_get_study_info(con,study_Id)
gs_lm_path <- study_info[['study_analysis_path']]
gs_lm_data <- suppressMessages(gs_lm_path %>%
	    gs_url() %>%
		gs_read() %>%
			select(ID,Active) %>% 
				filter(Active==1))
message("Active linear models:")
print(gs_lm_data)

lm_list <- gs_lm_data$ID
message("4. Getting latest results from database for these models")
print(site_VEC)
print(gs_lm_data$ID)
res<-eregr_meta_int_select_latest_lm_results(con,study_Id,site_VEC,lm_list)
saveRDS(res,file=paste(Results_Path,'/',"res_latest_lm.rds",sep=''))

if (n_param == "overall") {
	res <- res %>%
		filter(n_overall >= n_min) %>%
			select(-n_overall,-n_cont,-n_pat)
} else if (n_param == "cohd") {
	res <- res %>%
		filter((n_pat >= n_min) & (n_cont >= n_min) ) %>%
			select(-n_overall,-n_cont,-n_pat)
} else 
	stop("Incorrect name for minimum N parameter")

saveRDS(res,file=paste(Results_Path,'/',"res_latest_lm_filtered.rds",sep=''))

message("5. Running meta-analysis")
tic()
meta_sess_id_and_time <- eregr_int_get_unique_time_id()
meta_res <- map(ROI_VEC, function(x) {
                map(lm_list, safely(~eregr_meta_analysis_onemodel(con,study_Id,.,x,res,meta_sess_id_and_time,compare_models=TRUE)))
    })
toc()

saveRDS(meta_res,file=paste(Results_Path,'/',"res_meta.rds",sep=''))

message("6. Disconnecting from database")
is_disconnected <- eregr_disconnect(con)

message("Finished successfully")

#



#
#tic()
#meta_sess_id_and_time <- eregr_int_get_unique_time_id()
#map(ROI_list, function(x) {
#                map(lm_list, safely(~eregr_meta_analysis_onemodel(db_conn,study_Id,.,x,res,meta_sess_id_and_time)))
#    })
#toc()
