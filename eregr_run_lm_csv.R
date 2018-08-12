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
if(length(cmdargs) != 5) {
	dput(cmdargs)
	stop("incorrect number of parameters")
}

study_Id <- cmdargs[1]

site_Id <- cmdargs[2]


ROI_VEC_arg <- cmdargs[3]
ROI_VEC <- str_split(ROI_VEC_arg,",")[[1]]
names(ROI_VEC) <- ROI_VEC
#CURRENT_ROI <- cmdargs[3]

Results_Path <- cmdargs[4]

dbFile <- cmdargs[5]
subjectsExcl=NA

message("1. Session Info")
sessionInfo()


# Connect to db
message("2. Connecting to database")
source(file=dbFile)
con <- database_connect()

# check that study_exist
# check that site exist

message ("3. Running linear models")
tic()
df_list <- NULL 
lm_err <- try(df_list <- map(ROI_VEC, ~eregr_run_linear_models(con, study_Id, site_Id, . ,subjectsExcl)))
toc()
if(is.null(df_list))
	stop(lm_err)
message("4.1 Model computation finished. Brief summary:")

saveRDS(df_list,paste(Results_Path,"/dflist_",study_Id,"_",site_Id,".rds",sep=""))




err_byROI<-pmap(list(df_list_curroi=df_list,ROI_vec=ROI_VEC), function(df_list_curroi,ROI_vec) {
        errs <- map(df_list_curroi,"error")
        message("ROI: ", ROI_vec)
        for (i in seq_along(errs)) {
            if (is.null(errs[[i]]))
                message("\tSuccess on precheck: \t", names(errs)[[i]])
            else
                message("\tFailed on precheck: \t", names(errs)[[i]])
        }

        message("")
        results <- map(df_list_curroi,"result")
        res_null <- map_lgl(results, ~ is.null(.))
        results <- results[!res_null]
        results <- map(results,"lmres") # %>% transpose()
        errs_post <- pmap(list(res_by_model = results, res_names = names(results)),function(res_by_model,res_names) {
            output <- pmap(list(x=res_by_model,x_name = res_names), function (x,x_name) {
                    if (is.null(x[['error']])) {
                            message("\tSuccess on postcheck: \t",x_name,"\t", x[['metric']],"\t",x[['vertex']])
                            NULL
                    }
                    else {
                            message("\tFailed on postcheck: \t",x_name,"\t",x[['metric']],"\t",x[['vertex']])
                            x
                    }
            })    
            output
        })
        list(err_precheck = errs, err_postcheck = errs_post)
        })

names(err_byROI) <- ROI_VEC

message("5. Writing model results to database")

#print(str(df_list))
res <- tryCatch(map(ROI_VEC, ~ eregr_write_results_linear_models(con,study_Id,site_Id,.,df_list[[.]])),
	error = function (e) {
			message("Error: couldn't write results to database")
			stop(e)

		})
message("6. Disconnecting from database")
is_disconnected <- eregr_disconnect(con)

message("7. Saving errors, models and results to files")
saveRDS(err_byROI,paste(Results_Path,"/error_",study_Id,"_",site_Id,".rds",sep=""))

# if needed we can save all models for debug purposes
#saveRDS(df_list,paste("lm_",study_Id,"_",CURRENT_ROI,"_",site_Id,".rds",sep=""))

saveRDS(res,paste(Results_Path,"/res_",study_Id,"_",site_Id,".rds",sep=""))

message("Finished successfully")
