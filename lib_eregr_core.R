
library(dplyr)
library(RSQLite)
library(tidyr)
library(googlesheets)
library(lubridate)
library(DBI)
library(tictoc)
library(stringr)
library(purrr)
library(metafor)
library(forestplot)
#----auxiliary
var_as_num <- function (x,varname) {
    x[[varname]] <- as.numeric(x[[varname]])
    x
}
#---end auxiliary

partial.d<-function(t.val,df,n1,n2){
  d<-t.val*(n1+n2)/(sqrt(n1*n2)*sqrt(df))
  names(d)<-"effect size d"
  return(d)
}
CI1<-function(ES,se){
  ci<-c((ES-(1.96)*se),(ES+(1.96)*se))
  names(ci)<-c("95% CI lower","95% CI upper")
  return(ci)
}

se.d2<-function(d,n1,n2){
  se<-sqrt((n1+n2)/(n1*n2)+(d^2)/(2*(n1+n2-2)))
  names(se)<-"se for d"
  return(se)
}



#---error handling
condition <- function(subclass, par, message, call = sys.call(-1)) {
  structure(
    class = c(subclass, "error","condition"),
    list(message = message, call = call, parameters = par)    
  )
}



eregr_assert_model_precheck <- function(n_values) {
    assert_condition <- if (n_values[['fs_data_length']]==0)
	condition("model_precheck",n_values,"no record for model,ROI and metric in feature sets table. Running model is not possible") 
    else if (!is.na(n_values[['n_cont']]) & n_values[['n_cont']]<1) 
        condition("model_precheck",n_values,"no observations for control group")
    else if (!is.na(n_values[['n_pat']]) & n_values[['n_pat']]<1) 
        condition("model_precheck",n_values,"no observations for patient group")
    else if (!is.na(n_values[['cont_min']]) & !is.na(n_values[['n_cont']]) & n_values[['n_cont']]<n_values[['cont_min']])
        condition("model_precheck",n_values,"not enough observations for control group")
    else if (!is.na(n_values[['pat_min']]) & !is.na(n_values[['n_pat']]) & n_values[['n_pat']]<n_values[['pat_min']])
        condition("model_precheck",n_values,"not enough observations for patient group")
    else if (!is.na(n_values[['total_min']]) & !is.na(n_values[['n_total']]) & n_values[['n_total']]<n_values[['total_min']])
        condition("model_precheck",n_values,"not enough observations in total")
    else if (!is.na(n_values[['mainfactor_exists']]) & !is.na(n_values[['cohensD']]) & n_values[['mainfactor_exists']] & n_values[['cohensD']]) 
            if (!is.na(n_values[['levels']]) & n_values[['levels']]>2)
                condition("model_precheck",n_values,"main factor has more than 2 levels; could be a filtering problem")
            else if (!is.na(n_values[['levels']]) & n_values[['levels']]<2)
                condition("model_precheck",n_values,"main factor has less than 2 levels; could be a filtering problem")
            else if (!is.na(n_values[['n_cont']]) & !is.na(n_values[['n_pat']]) & 
                         n_values[['n_cont']]+n_values[['n_pat']]!=n_values[['n_total']])
                condition("model_precheck",n_values,"amount of controls and patients does not sum up to total; could be a filtering problem or problem with cont_val, pat_val")
            else
                NA
    else if( n_values[['mainfactor_exists']] & (!n_values[['cohensD']]) ) 
            if (n_values[['n_factors_in_model']]>0)
                condition("model_precheck",n_values,"variables for model contain factor. Partial correlations are not possible")
            else
                NA
	 
   else        
        NA
   if (mode(assert_condition)=="list")
        stop(assert_condition)
}
#---end error handling

eregr_write_results_linear_models <- function (db_conn, study_Id, site_Id, ROI,res_list) {
    eregr_int_write_results_one_linear_model <- function (db_conn, study_Id, site_Id, ROI, lm_name, lm_mainfactor, lm_res,
                                                            res_sessionID,
                                                          tbl_lmres = "lm_results",
                                                         tbl_lm_cohd_res = "lm_cohend_results",
                                                          tbl_lm_corr = "lm_corr_results",
                                                          tbl_lmvars = "lm_variables",
                                                          tbl_lm_res_keys="lm_results_keys",
							tbl_cont_summ="lm_summary_cont",
							tbl_fact_summ="lm_summary_fact",
							tbl_lm_demog_res="lm_demog_results") {

        
        # check if NULL
        if(is.null(lm_res)) stop("linear model ", lm_name, " finished with error")
        # process if not NULL
        # process main lm_results table
        res_demog <- lm_res[['precheck']]
        res_models <- lm_res[['lmres']]

	res_error <- map_lgl(res_models, function (res) {
			if (!is.null(res[['error']])){
#				message("ailed on postcheck: ", lm_name, "; metric: ", res[['metric']], "; vertex: ", res[['vertex']])
				TRUE
			}
			else
				FALSE
		})
	report_error <- res_models[res_error]



	res_models <- res_models[!res_error]
	


        lm_record <- eregr_get_lm_record(db_conn,lm_res[['sessionID']], lm_name)
        lm_mainfactor <- lm_record[['main_factor']]
        
        lm_ses_Id <- lm_res[['sessionID']]
        res_coef <- map(res_models, ~ coef(summary(.))[,1])
        res_sterr<- map(res_models, ~ coef(summary(.))[,2])
        res_pbeta <- map(res_models, ~ coef(summary(.))[,4])
        
        res_models_tp <- transpose (res_models)
        
        df_beta <- as.data.frame (do.call("rbind",res_coef))
        df_sterr <- as.data.frame (do.call("rbind",res_sterr))
        df_pbeta <- as.data.frame (do.call("rbind",res_pbeta))
        
        df_metr <- as.data.frame (do.call("rbind",res_models_tp[['metric']]))
        names(df_metr) <- 'metric'
#explicitly set metrics to be character
        df_metr[[1]] <- as.character(df_metr[[1]])
            
        df_vert <- as.data.frame (do.call("rbind",res_models_tp[['vertex']]))
        names(df_vert) <- 'vertex'
        
        df_beta <- cbind(df_beta,df_metr,df_vert)
        df_beta <- df_beta %>% gather (key="var",value="beta",-metric,-vertex)
        
        df_sterr <- cbind(df_sterr,df_metr,df_vert)
        df_sterr <- df_sterr %>% gather (key="var",value="sterr",-metric,-vertex)
        
        df_pbeta <- cbind(df_pbeta,df_metr,df_vert)
        df_pbeta <- df_pbeta %>% gather (key="var",value="p_beta",-metric,-vertex)
            
        df_res <- suppressMessages(inner_join(df_beta,df_sterr, by=c("metric","vertex","var")))
        df_res <- suppressMessages(inner_join(df_res, df_pbeta, by=c("metric","vertex","var")))
            
        #process lm_cohensd_results table
        
        #write main table lm_results    
        df_res[['lmID']] <- lm_name
        df_res[['sessionID']] <- lm_ses_Id
        df_res[['siteID']] <- site_Id
        df_res[['ROI']] <- ROI
        
        df_res <- df_res %>% select(lmID,sessionID,siteID,ROI,metric,vertex,var,beta,sterr,p_beta)
# extract and write keys into lm_results_keys table
        df_toreg_keys <- df_res %>% 
                        select(lmID,sessionID,metric,ROI) %>% 
                            distinct()
        metr_ul <- df_res %>%
                        select(metric) %>%	
				distinct() %>%
	                            unlist()
       l_write_keys <- dbWriteTable(db_conn, name = tbl_lm_res_keys, value = cbind(data.frame(res_keyID=NA),df_toreg_keys,data.frame(result_sessionID=res_sessionID)), append=TRUE,row.names=FALSE)

       if (l_write_keys == FALSE)
           stop(simpleError("could not write keys into lm_results_keys table"))
       query <- sprintf("SELECT * FROM %s WHERE result_sessionID='%s' AND lmID='%s' AND sessionID='%s' AND  ROI='%s'",tbl_lm_res_keys,res_sessionID, lm_name,lm_ses_Id,ROI);
       df_reg_keys <- dbGetQuery(db_conn, query)
 	
	l_write_cont_metr=TRUE
	l_write_fact_metr=TRUE
       #prepare summary statistics
       for (m_i in seq_along(metr_ul)) {
		cur_metr=try(which(res_models_tp[['metric']]==metr_ul[[m_i]])[[1]]) #first item in model results with current metrics
		if(class(cur_metr)=='try-error') next
		cur_model <- res_models_tp[['model']][[cur_metr]]
		lm_summary <- map(cur_model, ~summary(.))
		lm_sd<- map(cur_model,~tryCatch(sd(.), 
                		warning = function(x) NA,
		                error = function(x) NA)) # lm_sd is NA if variable is a factor

		lm_cont_summ_table <- pmap(list(elem=lm_summary[!is.na(lm_sd)],elem_sd=lm_sd[!is.na(lm_sd)],elem_name=names(lm_summary[!is.na(lm_sd)])), 
                                       function(elem,elem_sd,elem_name) {

                                        cont_val<-try(data.frame(mean=elem[['Mean']],median=elem[['Median']],Q1=elem[['1st Qu.']], 
                                                                Q3=elem[['3rd Qu.']],min=elem[['Min.']],max=elem[['Max.']],std=elem_sd,
                                                                var=elem_name))
                                        if(class(cont_val)=='try-error') {
                                            cont_val <- data.frame(mean=NA, median=NA, Q1=NA, Q3=NA, min=NA, max=NA,std=NA,var=elem_name)
                                        }
                                        cont_val

                })
		lm_cont_summ_table<-do.call("rbind",lm_cont_summ_table)
		lm_fact_summ_table <- pmap(list(elem=lm_summary[is.na(lm_sd)],elem_name=names(lm_summary[is.na(lm_sd)])), function(elem,elem_name) {

                                        fact_val<-tryCatch(data.frame(value=as.numeric(names(elem)),amount = elem,var=elem_name),
								warning=function(x) {message(str(elem))
											message(elem_name)
											}
								)
                                        if(class(fact_val)=='try-error') {
                                            fact_val <- data.frame(value=NA,amount=NA,var=elem_name)
                                        }
                                        fact_val

                })
		lm_fact_summ_table <- do.call("rbind", lm_fact_summ_table)
	
		if(length(lm_cont_summ_table)>0) {
			lm_cont_summ_table[['lmID']] <- lm_name
			lm_cont_summ_table[['sessionID']] <- lm_ses_Id
			lm_cont_summ_table[['ROI']] <- ROI
			lm_cont_summ_table[['metric']] <- metr_ul[[m_i]]
	
			lm_cont_summ_table <- lm_cont_summ_table %>%
                    			left_join(df_reg_keys, by=c("lmID","sessionID","ROI","metric")) %>%
						select(res_keyID,var,mean,median,Q1,Q3,min,max,std)
			l_write_cont_metr <- l_write_cont_metr & dbWriteTable(db_conn,name = tbl_cont_summ, value = lm_cont_summ_table, append=TRUE, row.names=FALSE)	
		}
		else
			l_write_cont_metr <- FALSE
		if(length(lm_fact_summ_table)>0) {
			lm_fact_summ_table[['lmID']] <- lm_name
			lm_fact_summ_table[['sessionID']] <- lm_ses_Id
			lm_fact_summ_table[['ROI']] <- ROI
			lm_fact_summ_table[['metric']] <- metr_ul[[m_i]]
			lm_fact_summ_table <- lm_fact_summ_table %>%
                    			left_join(df_reg_keys, by=c("lmID","sessionID","ROI","metric")) %>%
						select(res_keyID,var,value,amount)
		
			l_write_fact_metr <- l_write_fact_metr & dbWriteTable(db_conn, name = tbl_fact_summ, value = lm_fact_summ_table, append=TRUE, row.names = FALSE)
		}
		else
			l_write_fact_metr <- FALSE
	}
	
       df_res <- df_res %>% 
                    left_join(df_reg_keys, by=c("lmID","sessionID","ROI","metric")) %>%
                        select(res_keyID,vertex,var,beta,sterr,p_beta)

# join with new keys and remove redundant columns
            
#        df_demog <- data.frame (lmID = lm_name, sessionID = lm_ses_Id, studyID = study_Id, 
        l_write_lm_stats <- FALSE
        l_write_lm_results <- FALSE
	l_write_demog_stats <- FALSE
        #depends on type of mainfactor - write results table
        if (is.na(lm_mainfactor) | lm_mainfactor == "")
            l_write_lm_results <- dbWriteTable(db_conn, name = tbl_lmres, value = df_res, append=TRUE,row.names=FALSE)
	
        else if(str_detect(lm_mainfactor,"factor")==TRUE) {
            res_cohd <- map(res_models, function (x) {
                                            sum_x <- summary(x)
                                            res_tstat <- coef(sum_x)[,3]
                                            deg_fr <- sum_x[['df']][[2]]
                                            n_cont <- lm_res[['precheck']][['n_cont']]
                                            n_pat <- lm_res[['precheck']][['n_pat']]
					    n_overall <- lm_res[['precheck']][['n_total']]					    
                                            
					    lm_mf_idx <- which(str_detect(names(res_tstat),fixed(lm_mainfactor)))[[1]]
					   
                                            tstat_val <- res_tstat[[lm_mf_idx]]

                                            names_tstat<-names(res_tstat)
                                            var <- names_tstat[[lm_mf_idx]]
                                            coh_d <- partial.d (tstat_val,deg_fr,n_cont,n_pat)
                                            coh_se <- se.d2(coh_d,n_cont,n_pat)
                                            coh_bound <- CI1(coh_d,coh_se)
                                            coh_low_ci <- coh_bound[[1]]
                                            coh_high_ci <- coh_bound[[2]]
                                            coh_pval<- coef(sum_x)[lm_mf_idx,4]
                                            coh_res <- c(var,coh_d,coh_se,coh_low_ci,coh_high_ci,coh_pval,n_overall,n_cont,n_pat)
                                            names(coh_res) <- c("var","cohens_d","cohens_se","cohens_low_ci","cohens_high_ci","cohens_pval","n_overall","n_cont","n_pat")
                                            coh_res
                            })
            df_cohd <- as.data.frame (do.call("rbind",res_cohd))
            df_cohd <- cbind(df_cohd, df_metr,df_vert)
            df_cohd[['lmID']] <- lm_name
            df_cohd[['sessionID']] <- lm_ses_Id
            df_cohd[['siteID']] <- site_Id
            df_cohd[['ROI']] <- ROI
            
            df_cohd <- df_cohd %>%
                            left_join(df_reg_keys,by=c("lmID","sessionID","ROI","metric"))
	    df_demog <- df_cohd %>%
				select (res_keyID,n_overall,n_cont,n_pat)
	    df_cohd <- df_cohd %>%
                                select (res_keyID,vertex,var,cohens_d,cohens_se,cohens_low_ci,cohens_high_ci,cohens_pval)
            l_write_lm_stats <- dbWriteTable(db_conn, name = tbl_lm_cohd_res, value = df_cohd, append=TRUE,row.names=FALSE)

            l_write_demog_stats <- dbWriteTable(db_conn, name = tbl_lm_demog_res, value = df_demog, append=TRUE,row.names=FALSE)

            l_write_lm_results <- dbWriteTable(db_conn, name = tbl_lmres, value = df_res, append=TRUE,row.names=FALSE)
            
        }
        else {
           
	    # 1. check that mainfactor is the same as 1st variable
            res_pcor <- map (res_models, function (x) {
                      lm_mf_idx <- which(str_detect(names(x[['model']]),fixed(lm_mainfactor)))[[1]]
                      var <- names(x[['model']])[[lm_mf_idx]]
 		      n_cont <- 0
                      n_pat <- 0
		      n_overall <- lm_res[['precheck']][['n_total']]					    
                       
                      pcor <- try(ppcor::pcor.test(x$model[,1],x$model[,2],x$model[,3:ncol(x$model)]))
		      if (class(pcor)=="try-error"){
                      	corr_res <- c(var,-100,-100,-100,n_overall,n_cont,n_pat)
                    
                      	names(corr_res) <- c ("var","corr","corr_pval","corr_se","n_overall","n_cont","n_pat")
                      	return (corr_res)           
 
		      }

                      corr_val<-pcor[,1]
                      corr_pval<- pcor[,2]
     		                           

		     
                      corr_se <- summary(x)$coefficients[names(x[['model']])[[2]],2]
                      corr_res <- c(var,corr_val,corr_pval,corr_se,n_overall,n_cont,n_pat)
                    
                      names(corr_res) <- c ("var","corr","corr_pval","corr_se","n_overall","n_cont","n_pat")
                      corr_res              
            })
            df_pcor <- as.data.frame(do.call("rbind",res_pcor))
            df_pcor <- cbind (df_pcor,df_metr,df_vert)
            df_pcor[['lmID']] <- lm_name
            df_pcor[['sessionID']] <- lm_ses_Id
            df_pcor[['siteID']] <- site_Id
            df_pcor[['ROI']] <- ROI
            df_pcor <- df_pcor %>% 
                            left_join(df_reg_keys,by=c("lmID","sessionID","ROI","metric"))
            df_demog <- df_pcor %>%
				select (res_keyID,n_overall,n_cont,n_pat)
    	    df_pcor <- df_pcor %>%
	                        select(res_keyID,vertex,var,corr,corr_pval,corr_se)
            l_write_lm_stats <- dbWriteTable(db_conn,name = tbl_lm_corr,value = df_pcor, append=TRUE,row.names=FALSE)
	    
            l_write_demog_stats <- dbWriteTable(db_conn, name = tbl_lm_demog_res, value = df_demog, append=TRUE,row.names=FALSE)
            l_write_lm_results <- dbWriteTable(db_conn, name = tbl_lmres, value = df_res, append=TRUE,row.names=FALSE)
            
        }
        # think of the output
        list(lm_results=l_write_lm_results,lm_stats=l_write_lm_stats, lm_demog = l_write_demog_stats, lm_error = report_error)
    }
    res_list_tp <- transpose(res_list)
    model_res <- res_list_tp[['result']]
    lm_names <- names(res_list)
        
    res_id_and_time <- eregr_int_get_unique_time_id(db_conn)
    res_sessionID <- res_id_and_time[[2]]
       #temporarily removed 'safely'
    res<-pmap(list(lm_name=lm_names, lm_res = model_res), safely(eregr_int_write_results_one_linear_model), 
         db_conn = db_conn, study_Id = study_Id, site_Id = site_Id, ROI = ROI, res_sessionID = res_sessionID)
    names(res) <- lm_names
    res    
}


eregr_set_lm_results <- function (db_conn, study_Id, site_Id, ROI, lm_Id, lmres,tbl_lmres="lm_results") {
    write_line_lm_results <- function (study_Id,site_Id,ROI,lm_Id,lmres_item, metr_item, vert_item, header_lmres) {
        lm_summ <- summary(lmres_item)
        return (lm_summ)
        
    }
    header_lmres <- eregr_int_get_tbl_header(db_conn,tbl_lmres)

    lm_metrname <- map(lmres,"metric")
    lm_vertname <- map(lmres,"vertex")

    lmres_item <- lmres[[1]]
    metr_item <- lm_metrname[[1]]
    vert_item <- lm_vertname[[1]]

    write_line_lm_results(study_Id, site_Id, ROI, lm_Id, lmres_item, metr_item, vert_item, header_lmres)
}

eregr_set_lm_list_results <- function (db_conn, study_Id, site_Id, ROI, lmres_list) {
    
    
}

eregr_get_lm_vars <- function (db_conn,  session_Id, lm_Id = "", tbl_vars = "lm_variables",tbl_lm = "linear_model") {
    query_vars <- sprintf("SELECT * FROM %s WHERE sessionID='%s'",tbl_vars,session_Id)
    if(lm_Id != "") query_vars <- sprintf ("%s AND lmID='%s'",query_vars, lm_Id)
    res_vars <- dbGetQuery(db_conn, query_vars)
    query_lm <- sprintf ("SELECT lmID, lm_text FROM %s WHERE sessionID='%s'", 
                       tbl_lm, session_Id)
    if(lm_Id != "") query_vars <- sprintf ("%s AND lmID='%s'",query_lm, lm_Id)
    res_lm <- dbGetQuery(db_conn,query_lm)
    res <- suppressMessages(res_vars %>%
            left_join(res_lm))
    res %>%
        select (var, lm_text, lmID, modifier,is_global)
}
eregr_get_lm_mutate_vars <- function(db_conn, session_Id, lm_Id = "", tbl_mutate_vars = "lm_mutate") {
    query <- sprintf("SELECT * FROM %s WHERE sessionID='%s'",tbl_mutate_vars,session_Id)
    if (lm_Id != "") query <- sprintf("%s AND lmID='%s'", query, lm_Id)
    res <- dbGetQuery(db_conn, query)
    res %>%
        select (var, order, formula, lmID)    
}

eregr_get_lm_text <- function (db_conn, session_Id, lm_Id, tbl_linear_model = "linear_model") {
    query <- sprintf("SELECT lm_text FROM %s WHERE  lmID='%s' AND  sessionID='%s'",
                    tbl_linear_model, lm_Id, session_Id)
    df <- dbGetQuery(db_conn, query)
    return (df$lm_text)
    
}
eregr_get_lm_record <- function (db_conn, session_Id, lm_Id, tbl_linear_model = "linear_model") {
    query <- sprintf("SELECT * FROM %s WHERE lmID='%s' AND sessionID='%s'",
                    tbl_linear_model, lm_Id, session_Id)
    df <- dbGetQuery(db_conn, query)
    return (df)
    
}
    
eregr_get_lm_and_fs_ID <- function (db_conn, session_Id, tbl_linear_model = "linear_model") {
    query <- sprintf("SELECT lmID,fsID,session_fs_ID FROM %s WHERE sessionID='%s'",
                    tbl_linear_model, session_Id)
    df <- dbGetQuery(db_conn, query)
    res <-  transpose(list(df[['lmID']],df[['fsID']],df[['session_fs_ID']]))
    names(res) <- df [['lmID']]   
    return (res)
    
}

eregr_get_fs_covariates <- function (db_conn, study_Id, site_Id, fs_Id, session_fs_Id, tbl_fs_covars="feature_sets_covariates",tbl_sess_fs="session_fs", tbl_sites_in_study="sites_in_study") {
	query <- sprintf("SELECT FSC.session_fs_ID,FSC.fsID,subjID,cov_name,cov_value,metric,vertex FROM %s as FSC, %s as SFS, %s as SIS WHERE 
			FSC.session_fs_ID = SFS.session_fs_ID AND FSC.fsID = SFS.fsID AND SFS.study_site_ID= SIS.study_site_ID
			AND studyID = '%s' AND siteID='%s' AND FSC.session_fs_ID='%s' AND FSC.fsID='%s'
			",tbl_fs_covars,tbl_sess_fs,tbl_sites_in_study,study_Id,site_Id,session_fs_Id,fs_Id)
	res<-dbGetQuery(db_conn,query)
	res
}

eregr_run_linear_models <- function (db_conn, study_Id, site_Id, ROI , excludeSubj = NULL) {

    eval_lm_metrvert <- function (lm_text,metr_text,vert_text, df_data, var_list) {
        lm_fulltext <- paste("lm(res~",lm_text,",data = df_data)",sep="")
        res<- tryCatch( {
                lmres <- eval(parse(text=lm_fulltext))
        
                # post-check 1: amount of columns in model is equal to df_data:
                model_vars <- var_list %>% filter (modifier != "filter")
                n_vars <- nrow(model_vars) + 1
                ncol_data = ncol(df_data)
                ncol_model = ncol(lmres[['model']])
                names_data = names(df_data)
                names_model = names(lmres[['model']])
                nrow_data =nrow(df_data)
                nrow_model= nrow (lmres[['model']])
            
                assert_condition <- if( n_vars !=ncol(lmres[['model']]) )
                   condition("model_lmeval_postcheck",list(n_vars = n_vars, ncol_data = ncol_data, ncol_model = ncol_model, names_data = names_data, names_model = names_model, nrow_data = nrow_data, nrow_model = nrow_model ),"model takes less variables than exist in data. Maybe model is overdefined.")
                else if (nrow(df_data)!=nrow(lmres[['model']]))
                    condition("model_lmeval_postcheck",list(n_vars = n_vars, ncol_data = ncol_data, ncol_model = ncol_model, names_data = names_data, names_model = names_model, nrow_data = nrow_data, nrow_model = nrow_model),"not equal number of rows in input data and model. Presumably metrics (or covars) have NA values.")
                else 
                    NA
                if (mode(assert_condition)=="list")
                    stop(assert_condition)
                lmres[['error']]=NULL
                lmres
                },
                error = function(e) list('error'=e))
        res[['metric']] <- metr_text
        res[['vertex']] <- vert_text
        res
   }
   run_model <- function (lm_Id_elem,df_covar,metr_data) {
	#first check that there is metr data
	if(nrow(metr_data)==0)
		stop("No data in metrics table for study: ", study_Id, " site: ", site_Id,  " ROI: ", ROI)
	
	lm_Id <- lm_Id_elem[[1]]
	fs_Id <- lm_Id_elem[[2]]
	session_fs_Id <- lm_Id_elem[[3]]
	message("Running model ", lm_Id)
	global_fs_covars <- eregr_get_fs_covariates(db_conn,study_Id,site_Id,fs_Id,session_fs_Id)
#	saveRDS(global_fs_covars,'/ifs/loni/faculty/thompson/four_d/disaev/projects/mass_uv_regr_20/enigma_mdd/gl_fs_covar.rds')
	# augment metr_data with global covariates only if there are those covariates
	if (nrow(global_fs_covars)>0) {
		global_fs_covars <- global_fs_covars %>% spread(cov_name,cov_value)
		metr_data <- inner_join(metr_data, global_fs_covars, by=c("subjID","metric","vertex"))
	}
        var_list <- eregr_get_lm_vars(db_conn, session_Id, lm_Id)
        mutate_df <- eregr_get_lm_mutate_vars(db_conn, session_Id,lm_Id)
        mutate_df
        #mutate, select and cleanup
        CURRENT_ROI = ROI
        # transforming cov_value to numeric 
    
	df_covar[,'cov_value'] <- as.numeric (df_covar[,'cov_value'])
        # spreading
        df_covar <- df_covar %>%
                spread(key = cov_name, value = cov_value)

#	saveRDS(df_covar,'/ifs/loni/faculty/thompson/four_d/disaev/projects/mass_uv_regr_20/enigma_mdd/NESDA/res/df_covar_0.rds')
        #tocheck with Neda and Boris - filtering for complete cases
#        lmvars_plain <- intersect (var_list[['var']], names(df_covar))
#        df_covar <- df_covar[complete.cases(df_covar[,lmvars_plain]),]
        #end tocheck

#augment metr_data with new regressors
# SELECT FROM feature_sets_covars. Filters: studyId, siteId, fsID <- from lin.model, session_fs_ID <- from lin. model
	
        #mutating
        mutate_df <- mutate_df %>%
                           arrange(order)
	if(nrow(mutate_df)>0)
	        for (i in 1:nrow(mutate_df)) {
        	    df_covar[[mutate_df[i,'var']]] <- df_covar %>%
                                            with(eval(parse(text=mutate_df[i,'formula'])))
	        }
        #selecting
        try({   var_list_nonglob <- var_list %>% filter(is_global!=1)
		colNums <- match(var_list_nonglob[['var']],names(df_covar))
	     	df_covar <- df_covar %>%
			 select (colNums,subjID)
		var_list_notFilters <- var_list_nonglob %>% filter(modifier!='filter')
		colNumsNotFilt <-  match(var_list_notFilters[['var']],names(df_covar))
		df_covar <- df_covar %>%
		                   drop_na(colNumsNotFilt)                            
                            })

        lm_record <- eregr_get_lm_record(db_conn, session_Id, lm_Id)
        lm_text <- lm_record[['lm_text']]
        lm_mainfactor <- lm_record[['main_factor']]
 
#---begin refactor
#	else 
#		if (str_detect(lm_mainfactor,"factor")==TRUE) { #factor(...) or factor(...):var
#			str_match <- str_match_all(lm_mainfactor,"(factor)\\(([A-Za-z0-9]+)\\)")
#			factor <- as.data.frame(str_match[[1]])
#			if(nrow(factor)>1) {
#	        	    	assert_condition <- condition("model_precheck",factor,"main_factor consists of more than one factor variable")
#				stop(assert_condition)
#			}
#			else if (nrow(factor)==1) # it is 1 factor or factor:continuous interaction
#				#PotERR if we have only (factor) in main_factor text
#				main_factor <- as.character(factor[[3]])
#		}
#		else { #no "factor" in main_factor field
#				main_factor <- lm_mainfactor
#		}
#end refactor
        filter_text <- eregr_get_lm_filter(db_conn,  session_Id, lm_Id)

	        #filtering
        if(!is.na(filter_text) & filter_text != "")
            df_covar <- tryCatch(df_covar %>% 
                            filter_(filter_text),
				error=function(e) {
					message("error in filtering")
					message(filter_text)
					stop(e)
				})
        
#	saveRDS(df_covar,'/ifs/loni/faculty/thompson/four_d/disaev/projects/mass_uv_regr_20/enigma_mdd/NESDA/res/df_covar_1.rds')
  message("starting precheck")
	#pre-check should go here
        #building pre-check structure

 	main_factor <- NA
	mainfactor_exists <- FALSE
        cohensD <- FALSE
	factor <- NA
        if (!(is.na(lm_mainfactor) | lm_mainfactor=="" | is.null(lm_mainfactor)))
            mainfactor_exists <- TRUE

	if(mainfactor_exists) {
		# case 1 - lm_mainfactor has a word 'factor' in it
		if (str_detect(lm_mainfactor,"factor\\([A-Za-z0-9\\_]+\\)")==TRUE) {
			str_factor_match <- str_match_all(lm_mainfactor,"(factor)\\(([A-Za-z0-9\\_]+)\\)")[[1]]
			factor <- as.data.frame(str_factor_match)
			#check - no factor:factor integration is allowed
			if(nrow(factor)>1) {
	        	    	assert_condition <- condition("model_precheck",factor,"main_factor consists of more than one factor variable")
				stop(assert_condition)
			}
			else if (nrow(factor)==1) {# it is 1 factor or factor:continuous interaction
				main_factor <- as.character(factor[[3]])	
				cohensD=TRUE
			}
		}
		#case 2 - interaction of 2 continuous variables
		else if (str_detect(lm_mainfactor,"([A-Za-z0-9\\_]+)\\:([A-Za-z0-9\\_]+)")==TRUE) {
        	    	assert_condition <- condition("model_precheck",factor,"interaction of continuous variables in main_factor is not allowed")
			stop(assert_condition)
		}
		else
			main_factor <- lm_mainfactor			
	}

        factors_in_model <- FALSE
        n_factors_in_model <- 0
        pat_val <- NA
	cont_val <- NA
	n_pat <- NA
        n_cont <- NA
        levels <- NA
	
	
        if (cohensD) {
            levels <- n_distinct(df_covar[[main_factor]])
 	    cont_val <- if(is.na(lm_record[['cont_value']]) | is.null(lm_record[['cont_value']]) ) 0 else lm_record[['cont_value']] 
	    n_cont <- sum(df_covar[[main_factor]]==cont_val) 
      	    pat_val <- if(is.na(lm_record[['pat_value']]) | is.null(lm_record[['pat_value']]) ) 1 else lm_record[['pat_value']]
            n_pat <- sum(df_covar[[main_factor]]==pat_val) 
        }
        else {
            lm_vars <- eregr_get_lm_vars(db_conn,session_Id,lm_Id)
            lm_factors_num <- lm_vars %>%
                            filter(modifier=="factor") %>%
                                summarise(n=n())
            n_factors_in_model <- lm_factors_num[['n']]
        }
        cont_min <- lm_record[['cont_min']]
        pat_min <- lm_record[['pat_min']]
        total_min <- NA
            
        n_total <- nrow(df_covar)
        #precheck - either no main_factor or factor(main_factor) or partcor - then there should be no factors in model
        # or no mainfactor at all
        #assert pre-check conditions

	#get feature set for this:
	# - model
	# - metric
	# - ROI
	# filter the list of df_metr and df_vert and df_data according to feature set
	# RUN evak_lm_metrvert
	fs_data <- eregr_get_fs_for_lm_ROI(db_conn,fs_Id,session_fs_Id,ROI)
	fs_data_length <- nrow(fs_data)
        n_value <- data.frame(cont_min,pat_min,total_min,cont_val,pat_val,as.character(main_factor),n_pat,n_cont,n_total,levels,
                              mainfactor_exists, cohensD, n_factors_in_model,fs_data_length)
 
        eregr_assert_model_precheck(n_value)
        message("Passed precheck: ",lm_Id)
        #joining with metrics
        df_full <- suppressMessages(metr_data %>% 
                    inner_join(df_covar) %>%
                        rename(res = value))
        #splitting by metric and vertex
	
#	return(list(df_full,metr_data,df_covar))
        df_list <- df_full %>% split (list(df_full$metric,df_full$vertex))
        df_metrvert <- as.list(names(df_list)) %>%
                        map(strsplit,split="\\.")
        df_metrvert <- df_metrvert %>%
                            map(1) %>%
                                transpose()
        df_metr <- df_metrvert[[1]]
        df_vert <- as.integer(df_metrvert[[2]])
        
#_NEW_FEATURE_SET 
# here add computation of new covariates
#move it up to df_full
	metr_list <- fs_data %>% 
			select(metric) %>%
				distinct() %>%
					unlist()

	metr_list_bool <- df_metr %in% metr_list
	df_metr <- df_metr[metr_list_bool]
	df_vert <- df_vert[metr_list_bool]
	df_list <- df_list[metr_list_bool]
	
        res <- list()
        res[['lmres']] <- pmap(list(metr_text=df_metr,vert_text=df_vert, df_data = df_list),eval_lm_metrvert, lm_text = lm_text, var_list = var_list)
        res[['precheck']]<-n_value
        res[['main_factor']] <- main_factor
        res[['sessionID']] <- session_Id
        #post-check goes here
        model_vars <- var_list %>% filter (modifier != "filter")
        n_vars <- nrow(model_vars) + 1 # +1 means + 'res' column
        
            
        res
        
    }
    message("Current ROI: ", ROI)

    #0.read covariates. we may need it to drop variables with sd=0
    df_covar <- eregr_get_covariates(db_conn, study_Id, site_Id)
#    sd_df_covar <- sapply(df_covar,sd)
#    sd0_df_covar <- names(df_covar)[sd_df_covar==0 & !is.na(sd_df_covar)]
    #1. register lm_sessionID
    session_Id <- eregr_int_register_session_lm(db_conn,study_Id, site_Id)
    #2. read line from study_Id
    study_info <- eregr_get_study_info(db_conn,study_Id)
    gs_lm_path <- study_info[['study_analysis_path']]
    gs_lm_data <- suppressMessages(gs_lm_path %>%
                    gs_url() %>%
                        gs_read())
    #3. register models
    res_lm_reg <- pmap(list(lm_Id = gs_lm_data[['ID']]),safely(eregr_int_register_lm), db_conn = db_conn, session_Id = session_Id,study_Id = study_Id,
                          gs_data = gs_lm_data, gsheet_lm_path = gs_lm_path, df_covar=df_covar, use_uscores = TRUE)
    results <- transpose(res_lm_reg)[['result']]
    errs <- transpose(res_lm_reg)[['error']]
    if ( sum(map_lgl(errs,~ !is.null(.) )) > 0 ){
        assert_condition <- condition("lm_registration",errs,"Not all models could be registered. Error on the stage of linear model registration.")
#	saveRDS(res_lm_reg,file='lm_reg_error.rds')
        stop(assert_condition)
    }
    #4. process covariates

    if(!is.null(excludeSubj))
        df_covar <- df_covar %>% 
                        filter(!(subjID %in% unlist(excludeSubj)))
    #5. get metrics data
    metr_data <- switch(study_info[['study_data_format']],
                        'raw'=eregr_read_metrics_data_byROI(db_conn,study_Id, site_Id, ROI),
                        'csv'=eregr_get_metrics_data_byROI(db_conn,study_Id, site_Id, ROI))
#        if(study_info[['study_data_format']]=='raw')
#        metr_data <- eregr_read_metrics_data_byROI(db_conn,study_Id, site_Id, ROI)
    
    #6. map for all active linear models
    lm_Id_list <- eregr_get_lm_and_fs_ID(db_conn, session_Id) #Feature set will probably just go here - we need to select model according to feature set
    res <- map (lm_Id_list, safely(run_model),df_covar = df_covar,metr_data = metr_data)
    names(res) <- names(lm_Id_list)
    res
}



eregr_get_metrics_data_byROI <- function (db_conn, study_Id, site_Id, ROI, 
                                             tbl_metrics_data="metrics_data", tbl_sess_metrics="session_loadMetrics", tbl_sites_in_study="sites_in_study",tbl_study_metrics="study_metrics") {
    query <- sprintf("SELECT * 
                        FROM %s INNER JOIN %s USING (session_metr_ID) 
                        INNER JOIN %s USING (study_site_ID)
                        WHERE studyID='%s' AND siteID='%s' AND ROI = '%s' AND metric IN (
			SELECT metr_name AS metric FROM %s
			WHERE studyID='%s'
			)",
                    tbl_metrics_data, tbl_sess_metrics, tbl_sites_in_study,  study_Id,site_Id,ROI,tbl_study_metrics,study_Id)
    res <- dbGetQuery(db_conn,query)
    res
}

eregr_get_all_metrics_data <- function (db_conn, study_Id, site_Id,
                                             tbl_metrics_data="metrics_data", tbl_sess_metrics="session_loadMetrics", tbl_sites_in_study="sites_in_study",tbl_study_metrics="study_metrics") {
    query <- sprintf("SELECT * 
                        FROM %s INNER JOIN %s USING (session_metr_ID) 
                        INNER JOIN %s USING (study_site_ID)
                        WHERE studyID='%s' AND siteID='%s'  AND metric IN (
                        SELECT metr_name AS metric FROM %s
                        WHERE studyID='%s'
                        )",
                    tbl_metrics_data, tbl_sess_metrics, tbl_sites_in_study,  study_Id,site_Id,tbl_study_metrics,study_Id)
    res <- dbGetQuery(db_conn,query)
    res
}
                       
eregr_read_metrics_data_byROI <- function (db_conn, study_Id, site_Id, ROI, 
                                           tbl_metrics_file="metrics_file", tbl_metrics_data="metrics_data") {
    
    store_metr_shape_data <- function (f_path,subj_Id,metr_name,rec_data) {
        #read data from file 
        f_toRead <- file(f_path,"rb")
        cur_file_data<-readBin(f_toRead, numeric(), n= 999999, size=4)
        close(f_toRead)
        n_vert <- length(cur_file_data) #number of vertices
        if(n_vert<1) stop("number of vertices is equal to zero. file ",f_path," is corrupted")
        #write data to table
        rep_studyID <- rep (study_Id,n_vert)
        rep_siteID <- rep (site_Id, n_vert)
        rep_subjID <- rep (subj_Id, n_vert)
        rep_ROI <- rep (ROI, n_vert)
        rep_metric <- rep (metr_name, n_vert)
        vertex <- seq(from = 1, to = n_vert)
        value <- cur_file_data
        rec_data[1:n_vert,'studyID'] <- rep_studyID
        rec_data[1:n_vert,'siteID'] <- rep_siteID
        rec_data[1:n_vert,'subjID'] <- rep_subjID
        rec_data[1:n_vert,'ROI'] <- rep_ROI
        rec_data[1:n_vert,'metric'] <- rep_metric
        rec_data[1:n_vert,'vertex'] <- vertex
        rec_data[1:n_vert,'value'] <- value
        rec_data
#        l_write_shape_data <- dbWriteTable(db_conn, name = data_tblname, value = rec_data, append=TRUE,row.names=FALSE)
#        return (l_write_shape_data)
        
    }
    rec_data_empty <- eregr_int_get_tbl_header(db_conn,tbl_metrics_data)
    
    query <- sprintf("SELECT * FROM %s WHERE studyID='%s' AND siteID='%s' AND ROI = '%s'",
                    tbl_metrics_file,study_Id,site_Id,ROI)
    df_metrfiles <- dbGetQuery(db_conn,query)
    list_metrfiles <- split(df_metrfiles, seq(nrow(df_metrfiles)) )
    
    list_rec_data <- list_metrfiles %>% 
                        map (~store_metr_shape_data(.[['file_path']],.[['subjID']],.[['metric']],rec_data_empty) )
#    return(list_rec_data)
    df_rec_data <- plyr::rbind.fill(list_rec_data)
    df_rec_data <- df_rec_data %>% select (-row_ID)
    df_rec_data
#    res <- dbGetQuery(db_conn,query)
#    res
}                       
                       
                       
eregr_get_study_info <- function (db_conn, study_Id, study_tblname="studies") {
    study_tbl <- dbReadTable(db_conn, name = study_tblname)
    study <- study_tbl %>%
                filter(studyID == study_Id)
    return (study)
    
}

#eregr_register_models_within_session <- function (db_conn, study_Id,gsheet_lm_path,lm_Id='', session_Id='',  generate_session_Id = TRUE) {
#    if(generate_session_Id) session_Id <- eregr_int_register_session_lm(db_conn,study_Id)
    
#for debugging purposes - registering just one model
#    gs_data <- suppressMessages(gsheet_lm_path %>%
#                gs_url() %>%
#                    gs_read())
#    
#    if (lm_Id != '') {
#        eregr_int_register_lm(db_conn,lm_Id = lm_Id,session_Id = session_Id,study_Id = study_Id,
#                          gs_data = gs_data, gsheet_lm_path = gsheet_lm_path, use_uscores = TRUE)
#    }
#    else {
#        lm_list <- gs_data %>% 
#                        filter(Active == 1)
#        lm_list <- lm_list$ID
#        res <- rep(FALSE,length(lm_list))
#        for (i in seq_along(lm_list)) {
##            
#            res[i] <- eregr_int_register_lm(db_conn,lm_Id = lm_list[[i]], session_Id = session_Id, study_Id = study_Id,
#                                 gs_data = gs_data, gsheet_lm_path = gsheet_lm_path, use_uscores = TRUE)
#            if (res[i]==FALSE) 
#                print (paste("Model:", lm_list[i],"could not be loaded."))
#
#        }
#        return (lm_list[res])
#        return (sum(as.numeric(res))==length(res))
#    }
#}

eregr_get_lm_filter <- function (db_conn,  session_Id, lm_Id, tbl_filter = "lm_filter") {
    query <- sprintf("SELECT * FROM %s WHERE lmID='%s' and sessionID='%s'",
                    tbl_filter,  lm_Id, session_Id)
    df <- dbGetQuery(db_conn, query)
    return (df$filter_full)
    
}

#connect to database routine
eregr_connect <- function (db_path,db_class=RSQLite::SQLite(),user=NA,password=NA,dbname=NA,host=NA,port=NA) {
    #checking for existence of SQLite database:
    if (class(db_class)==class(RSQLite::SQLite())) {
	if (!file.exists(db_path))
		stop(paste("File ",db_path," does not exist. Cannot connect to a database",sep=""))
        return (dbConnect(db_class,db_path)) 
    }
    else if(class(db_class)==class(RMySQL::MySQL())) {
	return(dbConnect(db_class,user = user, password = password, dbname = dbname, host = host, port = port))
    }
}

eregr_disconnect <- function(db_conn) 
    return (dbDisconnect(db_conn))
 
eregr_src_sqlite <- function(db_path) {
    return (src_sqlite(db_path,create=FALSE))
    
}

eregr_int_get_tbl_header <- function (db_conn,tblname) {
    query <- paste("SELECT * FROM", tblname,"LIMIT 0")
    res <- dbSendQuery(db_conn,query)
    tbl <-dbFetch(res)
    dbClearResult(res)
    tbl
}
eregr_add_study_metric <- function(db_conn, study_Id, metr_name, tbl_study_metrics="study_metrics") {
	newmetr <- eregr_int_get_tbl_header(db_conn,tbl_study_metrics)
	newmetr[1,'studyID'] <- study_Id
	newmetr[1,'metr_name'] <- metr_name
	dbWriteTable(db_conn,name = tbl_study_metrics, value = newmetr, append=TRUE,row.names = FALSE)
}
eregr_remove_study_metric <- function(db_conn, study_Id, metr_name, tbl_study_metrics="study_metrics") {
	query <- sprintf("DELETE FROM %s WHERE studyID='%s' AND metr_name='%s'",
			tbl_study_metrics,study_Id,metr_name)
	dbExecute(db_conn,query)
}

eregr_register_study <- function (db_conn,study_Id, study_Name, gsheet_path, study_tblname="studies", study_metr_tblname="study_metrics") {

    #read sheet and select line with study_Id from it     
    gsheet <- suppressMessages(gsheet_path %>%
                gs_url())
    gs_model <- suppressMessages(gsheet %>%
                    gs_read() %>%
                        filter(ID==study_Id))
    #check if there's exactly one row with this ID
    if (nrow(gs_model) != 1) {
	err <- simpleError(paste("no rows in google sheet for the study: \t", study_Id,sep="")) 	
	stop(err)
    }    
    rec_study <- eregr_int_get_tbl_header(db_conn, study_tblname)
    #some renaming to fit the database naming convention
    rec_model <- gs_model %>%
                    rename(studyID=ID,study_analysis_path=AnalysisList_Path,study_demographics_path=DemographicsList_Path,study_data_format=Type,study_fs_path=FeatureSetList_Path)
    rec_model$study_gDoc_path=gsheet_path
    #selecting the traits to write in separae data table
    nTraits <- length(strsplit(gs_model$Trait,'[; ]+')[[1]])
    trait_names=paste("trait",c(1:nTraits),sep="_")
    rec_model <- rec_model %>%
                    separate(Trait,into=trait_names,sep='[; ]+')
    rec_model <- rec_model %>%
                    #gather traits:
                    gather(study_metrics_trait,metr_name,which(colnames(rec_model) %in% trait_names)) %>%
                        #drop study_metrics_trait
                        select(-study_metrics_trait)        
    rec_study_metrics <- rec_model %>% select (studyID,metr_name)
        
    #grouping other lines
    rec_model <- rec_model %>% select(-metr_name) %>% distinct()
    rec_study[1,'studyID'] <- rec_model[['studyID']]
    rec_study[1,'study_analysis_path'] <- rec_model[['study_analysis_path']]
    rec_study[1,'study_demographics_path'] <- rec_model[['study_demographics_path']]
    rec_study[1,'study_fs_path'] <- rec_model[['study_fs_path']]
    rec_study[1,'study_data_format'] <- rec_model[['study_data_format']]
    rec_study[1,'study_gDoc_path'] <- gsheet_path
    #writing to database
    l_studies_write <- dbWriteTable(db_conn,name = study_tblname, value = rec_study, append=TRUE, row.names = FALSE)
				
    if (!l_studies_write) {
	stop("could not write a study to 'studies' table: ",study_Id,"\n")
    }
    l_st_metrics_write <- tryCatch(dbWriteTable(db_conn,name = study_metr_tblname, value = rec_study_metrics, append=TRUE, row.names = FALSE),
				   error = function(e) FALSE)
    if (!l_st_metrics_write) {
	stop("could not write metrics information  into metrics table. Study: ", study_Id)
    }
    return (l_studies_write & l_st_metrics_write)
}
    
eregr_get_study_metrics <- function (db_conn, study_Id,study_metr_tblname="study_metrics") {
    tbl_metrics <- dbReadTable(db_conn, name = study_metr_tblname)
    study_metr <- tbl_metrics %>% 
                filter (studyID == study_Id)
    return (study_metr)
    
}
eregr_register_site <- function (db_conn,study_Id,site_Id,study_tblname="studies", site_tblname="sites_in_study") {
    #check if study exists in studies table
    tbl_studies <- dbReadTable(db_conn,name=study_tblname)
    study_exist <- tbl_studies %>%
                        filter(studyID==study_Id) %>%
                            summarise(n=n()) %>%
                                as.data.frame()
    if(study_exist$n!=1) stop("study ", study_Id, " does not exist in studies table")
    
    #append a new row to a sites table
    rec_site <- eregr_int_get_tbl_header(db_conn, site_tblname)
    rec_site[1,'studyID'] <- study_Id
    rec_site[1,'siteID'] <- site_Id
    rec_site[1,'study_site_ID'] <- NA #setting NA explicitly for auto-increment
    l_site_write <- tryCatch(dbWriteTable(db_conn,name=site_tblname,value=rec_site,append=TRUE,row.names=FALSE),
			     error = function(e) {message(e);FALSE})
    if (!l_site_write) stop ("couldn't site information to database: \t", site_Id)
    return (l_site_write)
}
    
eregr_int_get_unique_time_id <- function(db_conn=NA) {
    dt <- Sys.time()
    timestamp <- if (class(db_conn)=="MySQLConnection") format(as.POSIXlt(dt),'%Y-%m-%d %H:%M:%S') else dt
    unique_time_id <- paste(year(dt),sprintf("%02d",month(dt)),sprintf("%02d",day(dt)),hour(dt)*3600+minute(dt)*60+round(second(dt),digits=6),sep='')
    return (list(timestamp,unique_time_id))
}    
    
eregr_int_register_session_lm <- function(db_conn,study_Id,site_Id,session_tblname="session_lm_analysis") {
    
    id_and_time <- eregr_int_get_unique_time_id(db_conn)
    rec_session <- eregr_int_get_tbl_header(db_conn,session_tblname)
	
    rec_session[1,'session_analysis_ID'] <- id_and_time[[2]]
    rec_session[1,'studyID'] <- study_Id
    rec_session[1,'siteID'] <- site_Id
    rec_session[1,'timestamp'] <- id_and_time[[1]]
    
    l_session_write <- dbWriteTable(db_conn,name=session_tblname,value=rec_session,append=TRUE,row.names=FALSE)
    if(l_session_write != FALSE)
        return (id_and_time[[2]])
    else
        return (l_session_write)
}

eregr_int_register_lm <- function (db_conn, lm_Id, session_Id, study_Id, gs_data, gsheet_lm_path, df_covar,lm_tblname = "linear_model",vars_tblname ="lm_variables",
					tbl_session_fs="session_fs", tbl_sites_in_study="sites_in_study",tbl_sla="session_lm_analysis",  use_uscores = FALSE) {
    
    register_vars_lm <- function (vars_tblname ="lm_variables",inter_tblname="lm_interactions") {
        tbl_vars <- eregr_int_get_tbl_header(db_conn,vars_tblname)
        tbl_inter <- dbReadTable(db_conn,name=inter_tblname)

        lm_text <- rec_lm$lm_text
	var_global <- str_extract_all(raw_lm_text,"(\\{global\\:([\\w]+)\\})")
	var_global <- str_replace_all(var_global[[1]],"(\\{global\\:([\\w]+)\\})","\\2")

        #extracting variable names
        var_list <- strsplit(lm_text,"[:+ ]+")[[1]]
        unique_var_list <- unique(var_list)
        #extracting interactions:
        var_split <- strsplit(lm_text,"[+ ]+")[[1]]
        var1 <- c()
        var2 <- c()
        for (elem in var_split) {
	
            var_inter <- strsplit(elem, "[:]")[[1]]
            if (length(var_inter)==1) next
            else if (length(var_inter==2)) {
                var1 <- c(var1,var_inter[1])
                var2 <- c(var2,var_inter[2])
            }
            else {
                print("incorrect encoding of interaction - multiple colons (:) in one variable")
                return (FALSE)
            }
        }
        interactions <- data.frame(var1,var2,stringsAsFactors = FALSE)
        #extracting modifiers. for now we have only one modifier - 'factor'
        var_modifiers <- vector("character", length(unique_var_list))
        var_name <-  vector("character", length(unique_var_list))
	var_isglobal <- vector("character",length(unique_var_list))
        for (i in seq_along(unique_var_list)){
            has_factor <- (length(grep(x = unique_var_list[i],pattern = "factor\\(.+\\)"))==1)
            if (has_factor) {
                var_pure <- gsub(x =  unique_var_list[i],pattern = "factor\\((.+)\\)",replacement="\\1")
                var_modifiers[i] <- "factor"
                var_name[i] <- var_pure
	    	var_isglobal[i] <- var_pure %in% var_global 
            }
            else {
                var_modifiers[i] <- ""
                var_name[i] <-unique_var_list[i]
		
	    	var_isglobal[i] <- unique_var_list[i] %in% var_global
            }
        }   
        var_recs <- data.frame(var_name,var_modifiers,var_isglobal,stringsAsFactors = FALSE)
        #check if one variable exists with and without modifier - it should not be like that
        tbl_sum <- var_recs %>% 
                        group_by (var_name) %>%
                            summarise(n_repeats=n()) %>%
                                filter (n_repeats > 1)
        if (nrow(tbl_sum) > 0) {
            for (i in 1:nrow(tbl_sum)) {
                print (paste("Variable",tbl_sum$var_name,"is repeated with and without 'factor' modifier several times in linear model. Skipping the model", lm_Id))
            }
            return (FALSE)
        }
        #preparing to write variables to table lm_variables
        rec_lm_vars <- tbl_vars[0,]
        for (i in 1:nrow(var_recs)) {
            rec_lm_vars[i,'lmID'] <- lm_Id
            rec_lm_vars[i,'sessionID'] <- session_Id
            rec_lm_vars[i,'var'] <- var_recs$var_name[i]
            rec_lm_vars[i,'modifier'] <- var_recs$var_modifiers[i]
            rec_lm_vars[i,'is_global'] <- as.logical(var_recs$var_isglobal[i])
	}
        #writing variables to table 'lm_variables'
        rec_lm_vars <- rec_lm_vars %>% 
                            drop_na()
	l_var_write <-  dbWriteTable(db_conn,name=vars_tblname,value=rec_lm_vars,append=TRUE,row.names=FALSE)    

        #preparing to write interactions to table 'lm_interactions'
        rec_lm_inter <- tbl_inter[0,]
        if (nrow(interactions)>0)
        for (i in 1:nrow(interactions)) {
            rec_lm_inter[i,'lmID'] <- lm_Id
            rec_lm_inter[i,'sessionID'] <- session_Id
            rec_lm_inter[i,'var1'] <- interactions$var1[i]
            rec_lm_inter[i,'var2'] <- interactions$var2[i]
        }
        #writing interactions to table 'lm_interactions'
        l_inter_write <-  dbWriteTable(db_conn,name=inter_tblname,value=rec_lm_inter,append=TRUE,row.names=FALSE)    
	print (l_var_write & l_inter_write)
        return (if (l_inter_write & l_var_write) rec_lm_vars else rec_lm_vars[0,])
    }
            
    register_filter_lm <- function (filter_tblname = "lm_filter") {
        tbl_filters <- dbReadTable(db_conn,name=filter_tblname) 
        
        filter1 <- if (is.null(gs_lm$Filters_1) | is.na(gs_lm$Filters_1)) '' else gs_lm$Filters_1
        filter2 <- if (is.null(gs_lm$Filters_2) | is.na(gs_lm$Filters_2)) '' else gs_lm$Filters_2
        
        filter_full <- filter1
        if (filter2!='') 
            filter_full <- paste("(",filter1,") | (",filter2,")")        
        filter_full<- gsub(x = filter_full,pattern = "(__)",replacement = "")

        # add variables from filters to variables table
        vars_f1 <- unlist(regmatches(filter1,gregexpr("(__[a-zA-Z0-9_-]+__)",filter1)))
        vars_f2 <- unlist(regmatches(filter2,gregexpr("(__[a-zA-Z0-9_-]+__)",filter2)))
        vars_full <- c (vars_f1,vars_f2)
        vars_full <- gsub(x=vars_full, pattern="(__)",replacement="")
        vars_full <- unique(vars_full)

        query <- sprintf("SELECT * FROM %s WHERE lmID = '%s' AND sessionID = '%s'",
                         vars_tblname, lm_Id,session_Id)
        df_vars <- dbGetQuery(db_conn, query)
#        df_filter_vars <- data.frame(var = vars_full, modifier = 'filter')
        vars_toadd <- vars_full[ !(vars_full %in% df_vars$var)]
        rec_var <-df_vars[0,]
        for (i in seq_along(vars_toadd)) {
            rec_var[i,'var'] <- vars_toadd[[i]]
            rec_var[i,'modifier'] <- "filter"
            rec_var[i,'lmID'] <- lm_Id
            rec_var[i,'sessionID'] <- session_Id
        }
        dbWriteTable(db_conn, name = vars_tblname, value = rec_var, append=TRUE,row.names=FALSE)

        rec_filter <- tbl_filters[0,]
        rec_filter[1,]=vector("character",ncol(tbl_filters))
        rec_filter$lmID <- lm_Id
        rec_filter$sessionID <- session_Id
        rec_filter$filter_1 <- filter1
        rec_filter$filter_2 <- filter2
        rec_filter$filter_full <- filter_full
        l_filter_write <- dbWriteTable(db_conn,name = filter_tblname,value = rec_filter,append=TRUE,row.names=FALSE)
        return (l_filter_write)
    }
            
    register_mutate_lm <- function (mutate_tblname="lm_mutate") {
        tbl_mutate <- dbReadTable(db_conn,name = mutate_tblname)
        var_list <- rec_vars %>%
                        select (var)
        newRegressors_txt <- if (is.null(gs_lm$NewRegressors) | is.na(gs_lm$NewRegressors)) '' else gs_lm$NewRegressors
        newRegressors_txt <- str_replace_all(newRegressors_txt,"[\\n\\r\\s]+","")
	if (newRegressors_txt == '') return (TRUE)
        
        regr_split <- strsplit(newRegressors_txt,"[;]+")[[1]]
        rec_mutate <- tbl_mutate[0,]
        
        var_mutate <- vector ("character",0)
        formula_mutate <- vector ("character",0)
        for (i in seq_along(regr_split)) {
            if (use_uscores==TRUE) {
                regr_split[i]<- gsub(x = regr_split[i],pattern = "(__)",replacement = "")
            }
	    var1 = str_match(regr_split[i],"^(\\w+)\\=")
            rec_mutate[i,'var'] <- var1[1,2]

	    formula1 = str_match(regr_split[i],"\\=(.*)$")
            rec_mutate[i,'formula'] <- formula1[1,2]

            #check length of formula_split
#            if(length(formula_split)!=2) {
#                print (paste("incorrect formula for variable", formula_split[1],". exiting model."))
#                return (FALSE)
#            }
            if(grepl("eval",rec_mutate[i,'formula'], fixed=TRUE)==TRUE) {
#                    ROI <- "13"
#                    rec_mutate[i,'formula'] <- eval(rec_mutate[i,'formula']))
                    message("Modified formula for ", rec_mutate[i,'var'], ": ",rec_mutate[i,'formula'])
            }
            rec_mutate[i,'order'] <- i
            rec_mutate[i,'lmID'] <- lm_Id
            rec_mutate[i,'sessionID'] <- session_Id
        }
        l_mutate_write <- dbWriteTable(db_conn, name = mutate_tblname, value = rec_mutate, append=TRUE,row.names=FALSE)
    }

#    gsheet_lm <- gsheet_lm_path %>%
#                gs_url()
#    gs_lm <- gsheet_lm %>%
#                gs_read() %>%
    gs_lm <- gs_data %>%        
                    filter(ID==lm_Id)
    if (nrow(gs_lm) != 1) return(FALSE)
    if (gs_lm$Active !=1 ) {
        return(FALSE)
    }
    message ("Registering model: ",lm_Id)
    
    raw_lm_text<-gs_lm$LM
    tbl_lm <- dbReadTable(db_conn,name=lm_tblname)
    rec_lm <- tbl_lm[0,]
    rec_lm[1,] <- vector("character",ncol(tbl_lm))
    rec_lm$lmID <- lm_Id
    rec_lm$sessionID <- session_Id
    rec_lm$lm_gDoc_path <- gsheet_lm_path
    rec_lm$lm_name <- gs_lm$Name
    rec_lm$lm_text <- gsub("[[:space:]]", "", gs_lm$LM)
    rec_lm$lm_text <- str_replace_all(rec_lm$lm_text,"(\\{global\\:([\\w]+)\\})","\\2")
    rec_lm$main_factor <- gs_lm$MainFactor
    rec_lm$new_regressors <- gs_lm$NewRegressors
    rec_lm$cont_value <- gs_lm$ContValue
    rec_lm$pat_value <- gs_lm$PatValue
    rec_lm$cont_min <- gs_lm$ContMin
    rec_lm$pat_min <- gs_lm$PatMin
    rec_lm$comments <- gs_lm$Comments    
    
    #drop sd=0 for specified columns:
    sd0_targets <- str_extract_all(rec_lm$comments,"\\{drop_sd0\\:([\\w]+)\\}")[[1]]
    sd0_targets <- str_replace_all(sd0_targets,"(\\{drop_sd0\\:([\\w]+)\\})","\\2")
    for (i in seq_along(sd0_targets)){
	sd_covar <- df_covar %>%
			filter(cov_name==sd0_targets[[i]]) %>%
				summarise(sd=sd(cov_value))
#	sd_covar <- df_covar %>%
#			group_by(cov_name) %>%
#				summarise(sd=sd(cov_value)) %>%
#					filter(cov_name==sd0_targets[[i]])
	is_sd0=(!is.na(sd_covar[['sd']]) & sd_covar[['sd']]==0)
	if (is_sd0) {
#		message("Covar ", sd0_targets[[i]], " has sd=0")
		text1 <- str_replace_all(rec_lm$lm_text,sprintf("[\\+](%s|factor\\(%s\\))(\\+|$)",sd0_targets[[i]],sd0_targets[[i]]),"+")
		text2 <- str_replace_all(text1,sprintf("[\\+](%s|factor\\(%s\\))[\\:]",sd0_targets[[i]],sd0_targets[[i]]),"+")
		text3 <- str_replace_all(text2,sprintf("[\\:](%s|factor\\(%s\\))",sd0_targets[[i]],sd0_targets[[i]]),"")
		text4 <- str_replace_all(text3,"\\+$","")
	    	rec_lm$lm_text <- text4
	}
    }
    fsID <- gs_lm$FeatureSet
    query = sprintf("SELECT fsID,MAX(session_fs_ID) as max_sess_fs_ID,SIS.studyID,SIS.siteID
        FROM %s SF, %s SLA, %s SIS
        WHERE fsID='%s'
	AND session_analysis_ID='%s'
        AND SF.study_site_ID = SIS.study_site_ID
	AND SIS.studyID=SLA.studyID AND SIS.siteID=SLA.siteID
        GROUP BY fsID,studyID,siteID",tbl_session_fs,tbl_sla,tbl_sites_in_study,fsID,session_Id)
    sess_fs <- dbGetQuery(db_conn,query)
    if (nrow(sess_fs)!=1) 
	stop("couldn't retrieve feature_set session!")
    rec_lm$fsID <- fsID
    rec_lm$session_fs_ID <- sess_fs$max_sess_fs_ID

    l_lm_write <- dbWriteTable(db_conn,name=lm_tblname,value=rec_lm,append=TRUE,row.names=FALSE)    
    
    rec_vars <- register_vars_lm()
    if (nrow(rec_vars) == 0) stop(paste("Could not register variables for model ", lm_Id,sep=""))
    rec_filter <- register_filter_lm()
    if (rec_filter == FALSE) stop(paste("Could not register filters for model ", lm_Id, sep=""))
    
    rec_mutate <- register_mutate_lm()
    if (rec_mutate == FALSE) stop(paste("Could not register new regressors for model ", lm_Id, sep=""))

    return (l_lm_write)
}

eregr_read_shape_all_subjects <- function (db_conn, data_dir, subj_list, roi_list, metrics_list,
                                    study_Id, site_Id, exclude_subj_file="", exclude_filter = 2,
                                                template = "<metric>_<ROI>.raw", rewrite = TRUE, metrfile_tblname = "metrics_file", exclude_subj_tblname = "subj_excluded") {
    
    
    erase_site_data <- function() {
        query <- sprintf("DELETE FROM %s WHERE studyID='%s' AND siteID='%s'",
                        metrfile_tblname,study_Id,site_Id)
        dbExecute(db_conn,query)
        query <- sprintf("DELETE FROM %s WHERE studyID='%s' AND siteID='%s'",
                        exclude_subj_tblname,study_Id,site_Id)
        dbExecute(db_conn,query)

    }
    if (rewrite)
        erase_site_data()

    exclude_subj<-data.frame()
    if (exclude_subj_file != "") {
        exclude_subj <- read.csv(exclude_subj_file,stringsAsFactors = FALSE)
        roi_cols <- which(str_detect(names(exclude_subj),"ROI"))
        exclude_subj <- exclude_subj %>% 
                            gather (key = "ROI", value = "QC", roi_cols)
        exclude_subj <- exclude_subj %>% 
                            rename(subjID = SubjID) %>%
                                select (subjID, ROI, QC)
        exclude_subj[['ROI']]<-str_replace(exclude_subj[['ROI']],"ROI","")
    }
    ses_metr_ID<-eregr_int_register_session_loadmetr(db_conn,study_Id,site_Id)
    res<-map(subj_list,~eregr_read_shape_one_subject(db_conn,data_dir,roi_list,metrics_list,study_Id,.,site_Id,ses_metr_ID, rewrite = rewrite) )
    #remove all empty results (those that finished correctly)
    res <- res[sapply(res,length)>0]
    if(length(res)==0) return(TRUE)
    #if there is subjects to exclude - continue and write it to database
    res <- do.call("rbind",res)
    res <- res %>% 
               select(subjID,roi_list) %>%
                    distinct()
    #read header of exclude_subj table:
    if(nrow(res)>0)
        for (i in 1:nrow(res)) {
            exclude_subj[['QC']] <- exclude_subj %>% with (replace (QC, (ROI == res[i,'roi_list'] & subjID == res[i,'subjID']),0)  )

        }
    exclude_subj[['studyID']] <- study_Id
    exclude_subj[['siteID']] <- site_Id
    exclude_subj[['session_metr_ID']] <- ses_metr_ID
    exclude_subj <- exclude_subj %>% select(siteID,studyID,subjID,ROI,QC,session_metr_ID)
    l_read_all_subj <- dbWriteTable(db_conn,name = exclude_subj_tblname, value = exclude_subj, append=TRUE,row.names=FALSE)
    l_read_all_subj
}

# this section will be devoted to metrics loading

eregr_int_register_session_loadmetr <- function(db_conn,study_Id,site_Id,session_tblname="session_loadMetrics",tbl_sites_in_study="sites_in_study") {
    
    query=sprintf("SELECT * FROM %s WHERE siteID='%s' AND studyID='%s'",
                   tbl_sites_in_study,site_Id,study_Id);
    rec_sites_in_study <- dbGetQuery(db_conn,query)
    
    id_and_time <- eregr_int_get_unique_time_id(db_conn)
    rec_session_metr <- eregr_int_get_tbl_header(db_conn,session_tblname)
    rec_session_metr[1,'session_metr_ID'] <- id_and_time[[2]]
    rec_session_metr[1,'timestamp'] <-id_and_time[[1]]
    if(nrow(rec_sites_in_study)==0)
        stop (simpleError("Site is not registered for the study"))
    else if (nrow(rec_sites_in_study)>1) 
        stop(simpleError("multiple rows for site and study. Database inconsistency"))
    rec_session_metr[1,'study_site_ID'] <- rec_sites_in_study[['study_site_ID']]
    l_session_write <- dbWriteTable(db_conn,name=session_tblname,value=rec_session_metr,append=TRUE,row.names=FALSE)
    if(l_session_write != FALSE)
        return (id_and_time[[2]])
    else
        return (l_session_write)
}

eregr_read_shape_one_subject <- function(db_conn, data_dir, roi_list, metrics_list, 
                                    study_Id, subj_Id,site_Id,ses_metr_Id,
                                    template = "<metric>_<ROI>.raw", rewrite = TRUE, metrfile_tblname = "metrics_file") {
    
    add_fname <- function (metr_val,roi_val,subj_Id,data_dir, template,rec_metr_file) {
        
            f_name <- gsub(template, pattern = "<metric>",replacement = metr_val)
            f_name <- gsub(f_name, pattern = "<ROI>", replacement = roi_val)
            f_path <- paste(data_dir,"/",subj_Id,"/",f_name,sep="")
            if (!file.exists(f_path)) {
                message("File",f_path,"does not exist. Data can not be read.")
                stop("File does not exist")
            }
            
        
            f_data <- file.info(f_path)
            f_lastedit <- f_data$mtime
            
            if(f_data$size==0) {
                return (FALSE)
                #stop("File ",f_path, " has size 0. Data can not be read")
            }
        
            
            rec_metr_file[1,'studyID'] <- study_Id
            rec_metr_file[1,'siteID'] <- site_Id
            rec_metr_file[1,'subjID'] <- subj_Id
            rec_metr_file[1,'metric'] <- metr_val
            rec_metr_file[1,'file_path'] <- f_path
            rec_metr_file[1,'ROI'] <- roi_val
            rec_metr_file[1,'session_metr_ID'] <- ses_metr_Id
            rec_metr_file[1,'file_lastedit'] <- f_lastedit
            rec_metr_file
        
    }

    rec_metr_file_empty <- eregr_int_get_tbl_header(db_conn,metrfile_tblname)

    df_metr <- data.frame(metrics_list, stringsAsFactors = FALSE)
    df_roi <- data.frame(roi_list, stringsAsFactors = FALSE)
    df_metroi <- merge(x=df_metr,y=df_roi,by=NULL)

    df_metroi_input <- split(df_metroi, seq(nrow(df_metroi)) )
    df_metroi <- map(df_metroi_input,~add_fname(.[[1]],.[[2]],subj_Id,data_dir,template,rec_metr_file_empty))
    exclude_metroi<-list()
    if (sum(unlist(df_metroi) == FALSE)>0) {
        res_exclude <- sapply(df_metroi,function(x) {if (is.logical(x)) 
                                !x
                            else 
                                FALSE })
        exclude_metroi <- df_metroi_input[res_exclude]
        df_metroi <- df_metroi[!res_exclude]
    }
    rec_metr_file <- do.call("rbind",df_metroi)
    l_write_metr_file <- dbWriteTable(db_conn,name = metrfile_tblname, value = rec_metr_file, row.names = FALSE,append = TRUE )
    #rbind exclude and write to exclude_subj folder
    rec_exclude_metroi <- do.call("rbind",exclude_metroi)
    if(!is.null(rec_exclude_metroi)) rec_exclude_metroi$subjID <- subj_Id
    rec_exclude_metroi                                                   
}

eregr_int_exists_covars <- function (db_conn, study_Id, site_Id, 
                                       tbl_covariates = "covariates_general", tbl_sites_in_study="sites_in_study",
                                        tbl_session_cov ="session_covariates") {
    query <- sprintf("SELECT *
                        FROM %s LEFT JOIN %s USING (study_site_ID)
                        INNER JOIN %s USING (session_covar_ID)
                        WHERE siteID='%s' AND studyID='%s'",
                        tbl_sites_in_study,tbl_session_cov,tbl_covariates, site_Id,study_Id)
    
    tbl_metrdata_record <- dbGetQuery(db_conn,query)
    return (nrow(tbl_metrdata_record))
}
eregr_int_erase_covars <- function (db_conn, study_Id, site_Id, 
                                       tbl_covariates = "covariates_general",tbl_session_cov="session_covariates",
                                       tbl_sites_in_study="sites_in_study") {
    
    query <- sprintf("DELETE FROM %s WHERE session_covar_ID IN 
                        (SELECT session_covar_ID
                            FROM %s INNER JOIN %s USING (study_site_ID)
                            WHERE siteID='%s' AND studyID='%s')",tbl_covariates,tbl_session_cov, tbl_sites_in_study,site_Id, study_Id)

    dbExecute(db_conn,query)
}
eregr_int_register_covar_session <- function (db_conn, study_Id, site_Id, cov_file_path, 
                                       tbl_sess_covariates = "session_covariates", tbl_sites_in_study="sites_in_study") {
        
        rec_covar_session <- eregr_int_get_tbl_header(db_conn, tblname = tbl_sess_covariates)
        time_and_id <- eregr_int_get_unique_time_id(db_conn)
        rec_covar_session[1,'session_covar_ID'] <- time_and_id[[2]]
        rec_covar_session[1,'timestamp'] <- time_and_id[[1]]
        if (cov_file_path=="update") {
            cov_file_path = as.character(as.POSIXct(time_and_id[[1]],origin="1970-01-01"))
            print(cov_file_path)
            rec_covar_session[1,'file_lastedit'] <- time_and_id[[1]]
        }
        else
            rec_covar_session[1,'file_lastedit'] <- file.info(cov_file_path)['mtime']
        rec_covar_session[1,'file_cov_path'] <- cov_file_path
                
        query = sprintf("SELECT * FROM %s WHERE siteID='%s' AND studyID='%s'",tbl_sites_in_study,
                            site_Id, study_Id);
        rec_sites_in_study <- dbGetQuery(db_conn,query);
        if(nrow(rec_sites_in_study)==0)
            stop (simpleError("Site is not registered for the study"))
        else if (nrow(rec_sites_in_study)>1) 
            stop(simpleError("multiple rows for site and study. Database inconsistency"))
        
        rec_covar_session[1,'study_site_ID'] <- rec_sites_in_study[['study_site_ID']] 
        l_write_covar_session <- dbWriteTable(db_conn, name = tbl_sess_covariates, 
                                              value = rec_covar_session, append=TRUE,row.names=FALSE)
        
        if (l_write_covar_session) 
            return (time_and_id[[2]]) 
        else {
            print("Could not register covariates session in database")
            return (FALSE)
        }
    }

    eregr_int_read_covariates <- function (cov_file_path) {
        cov_data <- read.csv(cov_file_path,dec = ".", sep=",", stringsAsFactors = FALSE)
        for (i in seq_along(cov_data)) {
            if (colnames(cov_data)[[i]]=="SubjID") next
            else 
                cov_data[,i] <- tryCatch( {
                        as.numeric(cov_data[,i])
                    },
                    warning = function(w) {
                        message ("Warning: possible problems with conversion of data to numeric. Please check your data.")
                        as.numeric(cov_data[,i])
#                        rep(c("NULL"),nrow(cov_data))
                    },
                    error = function(e) {
                        message ("Error: possible problems with conversion of data to numeric. Please check your data.")
                        message(e)
                        rep(c("NULL"),nrow(cov_data))
                },
                    finally = {
                    }
                )
                if ("NULL" %in% cov_data[,i]) {
                    message ("Exiting function because of possible problems with conversion.")
                    return (cov_data[0,])
                }

        }
        
        return (cov_data)
    }

            
eregr_register_covariates <- function (db_conn, study_Id, site_Id, cov_file_path, 
                                       tbl_covariates = "covariates_general", tbl_sess_covariates = "session_covariates") {

    
    
    cov <- eregr_int_read_covariates(cov_file_path)
    if (nrow(cov)<1) stop("Could not read covariates")

    if(eregr_int_exists_covars(db_conn, study_Id, site_Id, tbl_covariates)>0) 
        eregr_int_erase_covars (db_conn, study_Id, site_Id, tbl_covariates)

    cov_sess_ID <- eregr_int_register_covar_session (db_conn, study_Id, site_Id, cov_file_path, 
                                       tbl_sess_covariates)
    if (cov_sess_ID == FALSE) stop ("Could not register session for covariates")
     
    #edit - it should be data frame with zero rows.
    
    rec_tbl_cov <- eregr_int_get_tbl_header(db_conn, tblname = tbl_covariates)
    cov_gathered <- cov %>%
                        gather(key = "cov_name", value = "cov_value", -SubjID)
    cov_gathered <- cov_gathered %>%
                        rename(subjID = SubjID)
    rec_tbl_cov[1:nrow(cov_gathered),'subjID'] <- cov_gathered$subjID
#    rec_tbl_cov[1:nrow(cov_gathered),'studyID'] <- rep(study_Id,nrow(cov_gathered))
#    rec_tbl_cov[1:nrow(cov_gathered),'siteID'] <- rep(site_Id,nrow(cov_gathered))
    rec_tbl_cov[1:nrow(cov_gathered),'cov_name'] <-cov_gathered$cov_name
    rec_tbl_cov[1:nrow(cov_gathered),'cov_value'] <-cov_gathered$cov_value
    rec_tbl_cov[1:nrow(cov_gathered),'session_covar_ID'] <- rep(cov_sess_ID,nrow(cov_gathered))

    l_tbl_cov_write <- dbWriteTable(db_conn, name = tbl_covariates, value = rec_tbl_cov, append=TRUE,row.names=FALSE)
    if (l_tbl_cov_write) return (rec_tbl_cov) else stop("Could not write covariates data to the database table")
}

eregr_int_exists_metr <- function (db_conn, study_Id, site_Id, metric_name,
                                       tbl_metrics = "metrics_data", tbl_sites_in_study="sites_in_study",
                                        tbl_session_metr ="session_loadMetrics") {
    query <- sprintf("SELECT *
                        FROM %s LEFT JOIN %s USING (study_site_ID)
                        INNER JOIN %s USING (session_metr_ID)
                        WHERE siteID='%s' AND studyID='%s' AND metric='%s'",
                        tbl_sites_in_study,tbl_session_metr,tbl_metrics, site_Id,study_Id, metric_name)
    tbl_metrdata_record <- dbGetQuery(db_conn,query)
    return (nrow(tbl_metrdata_record))
}
eregr_int_erase_metr <- function (db_conn, study_Id, site_Id,metric_name, 
                                       tbl_metrics = "metrics_data", tbl_sites_in_study="sites_in_study",
                                        tbl_session_metr ="session_loadMetrics") {
    
    query <- sprintf("DELETE FROM %s WHERE session_metr_ID IN 
                        (SELECT session_metr_ID
                            FROM %s INNER JOIN %s USING (study_site_ID)
                            WHERE siteID='%s' AND studyID='%s' AND metric='%s')",
			tbl_metrics,tbl_session_metr, tbl_sites_in_study,site_Id, study_Id,metric_name)

    dbExecute(db_conn,query)
}


eregr_register_metrics_csv <- function (db_conn, study_Id, site_Id, metr_file_path, metric_name,
                                        tbl_metrics = "metrics_data", 
                                        tbl_sess_metrics="session_loadMetrics", 
                                        tbl_sites_in_study="sites_in_study") {
    metr <- eregr_int_read_covariates(metr_file_path)
    if (nrow(metr)<1) stop("Could not read metrics")
   if(eregr_int_exists_metr(db_conn, study_Id, site_Id, metric_name)>0) 
        eregr_int_erase_metr (db_conn, study_Id, site_Id,metric_name)
    metr_sess_ID <- eregr_int_register_session_loadmetr(db_conn,study_Id,site_Id)
    if (metr_sess_ID == FALSE) stop ("Could not register session for metrics data")

    rec_tbl_metr <- eregr_int_get_tbl_header(db_conn, tblname = tbl_metrics)
    metr_gathered <- metr %>% 
                        gather(key = "ROI", value="value", -SubjID)
    metr_gathered <- metr_gathered %>%
                        rename(subjID = SubjID)

    rec_tbl_metr[1:nrow(metr_gathered),'session_metr_ID'] <- rep(metr_sess_ID,nrow(metr_gathered))
    rec_tbl_metr[1:nrow(metr_gathered),'subjID'] <- metr_gathered$subjID
    rec_tbl_metr[1:nrow(metr_gathered),'metric'] <- rep(metric_name,nrow(metr_gathered))
    rec_tbl_metr[1:nrow(metr_gathered),'ROI'] <-metr_gathered$ROI
    rec_tbl_metr[1:nrow(metr_gathered),'vertex'] <- rep(1,nrow(metr_gathered))
    rec_tbl_metr[1:nrow(metr_gathered),'value'] <- metr_gathered$value
    rec_tbl_metr[1:nrow(metr_gathered),'row_ID'] <- rep(NA,nrow(metr_gathered))
    l_tbl_cov_write <- dbWriteTable(db_conn, name = tbl_metrics, value = rec_tbl_metr, append=TRUE,row.names=FALSE)
    if (l_tbl_cov_write) return (rec_tbl_metr) else stop("Could not write covariates data to the database table")

}        
        
        
eregr_get_subjlist_from_covars <- function (db_conn, study_Id, site_Id, tbl_covariates = "covariates_general",tbl_sess_covariates = "session_covariates",tbl_sites_in_study="sites_in_study") {
    #select the whole subtable of covariates for study/site from the covariates table. 
    #It's not the most efficient way.
    query <- sprintf("SELECT studyID,siteID,subjID,cov_name,cov_value,session_covar_ID
                            FROM %s INNER JOIN %s USING (session_covar_ID)
                            INNER JOIN %s USING (study_site_ID)
                            WHERE studyID='%s' AND siteID='%s'",tbl_covariates,tbl_sess_covariates, tbl_sites_in_study, study_Id,site_Id)
    tbl_covars <- dbGetQuery(db_conn, statement = query)
    df_subj <- tbl_covars %>%
                    select (subjID) %>% 
                        distinct()
    return (df_subj)
}

eregr_get_covariates <- function (db_conn, study_Id, site_Id, tbl_covariates = "covariates_general",tbl_sess_covariates = "session_covariates",tbl_sites_in_study="sites_in_study") {
    query <- sprintf("SELECT studyID,siteID,subjID,cov_name,cov_value,session_covar_ID  
                            FROM %s INNER JOIN %s USING (session_covar_ID)
                            INNER JOIN %s USING (study_site_ID)
                            WHERE studyID='%s' AND siteID='%s'",tbl_covariates,tbl_sess_covariates, tbl_sites_in_study, study_Id,site_Id)
    tbl_covars <- dbGetQuery(db_conn, statement = query)    
    return (tbl_covars)
}

eregr_export_covariates <- function (db_conn, study_Id, site_Id, output_path,tbl_covariates = "covariates_general") {
    tbl_covars <- eregr_get_covariates(db_conn, study_Id, site_Id, tbl_covariates)
    tbl_forexport <- tbl_covars %>%
                        select(-studyID,-siteID, -session_covar_ID) %>%
                            spread (cov_name,cov_value)
    session <- tbl_covars %>%
                    select (session_covar_ID) %>%
                        distinct()
    if (nrow(session)!=1) {
        print ("Data is corrupted, several sessions per one site and study")
        return (FALSE)
    }
    fname <- paste("covars",study_Id, site_Id,session$session_covar_ID,sep="_")
    fname <- paste(output_path,"/",fname,".csv",sep="")
    write.csv(file = fname, x = tbl_forexport)
    return (tbl_forexport)
    
}

#----- META-ANALYSIS FUNCTIONALITY
        eregr_int_select_active_models <- function(db_conn, study_Id) {
    study_info <- eregr_get_study_info(db_conn,study_Id)
    gs_lm_path <- study_info[['study_analysis_path']]
    gs_lm_data <- suppressMessages(gs_lm_path %>%
                    gs_url() %>%
                        gs_read())
    gs_lm_active <- gs_lm_data %>% 
                        filter(Active==1)
    gs_lm_id <- gs_lm_active[['ID']]
    gs_lm_id
}

eregr_meta_compare_mod_sess_pair <- function(db_conn, lm_Id, ROI, site_Id_1,site_Id_2,session_Id_1, session_Id_2, 
                                             tbl_linear_model="linear_model", tbl_lm_vars="lm_variables", 
                                            tbl_lm_filter="lm_filter", tbl_lm_inter="lm_interactions",
                                            tbl_lm_mutate="lm_mutate") {
    res<- list()
    result<-list()
#1. compare main linear_model tables
    query_1 <- sprintf("SELECT * FROM %s WHERE lmID='%s' AND sessionID='%s'",tbl_linear_model,lm_Id,session_Id_1)
    query_2 <- sprintf("SELECT * FROM %s WHERE lmID='%s' AND sessionID='%s'",tbl_linear_model,lm_Id,session_Id_2)
    
    df_lm_1 <- dbGetQuery(db_conn, query_1)
    df_lm_2 <- dbGetQuery(db_conn, query_2)

    df_lm_1_exact <- df_lm_1 %>% select(lmID,lm_gDoc_path,lm_name,cont_value,pat_value,cont_min,pat_min)
    df_lm_2_exact <- df_lm_2 %>% select(lmID,lm_gDoc_path,lm_name,cont_value,pat_value,cont_min,pat_min)
    l_main_lm_compare <- setequal(df_lm_1_exact,df_lm_2_exact)
    res[['lm']] <- if(l_main_lm_compare==TRUE) TRUE else list(df_lm_1_exact,df_lm_2_exact)
#    message("linear_model compare: ", l_main_lm_compare)
#    if(!l_main_lm_compare) return(FALSE)  #in next versions return conditions
#2. compare lm_variables
    query_1 <- sprintf("SELECT * FROM %s WHERE lmID='%s' AND sessionID='%s'",tbl_lm_vars,lm_Id,session_Id_1)
    query_2 <- sprintf("SELECT * FROM %s WHERE lmID='%s' AND sessionID='%s'",tbl_lm_vars,lm_Id,session_Id_2)
    
    df_var_1 <- dbGetQuery(db_conn, query_1)
    df_var_2 <- dbGetQuery(db_conn, query_2)
        
    df_var_1_exact <- df_var_1 %>% 
                        select(lmID,var,modifier) %>% 
                            arrange(lmID,var,modifier)
    df_var_2_exact <- df_var_2 %>% 
                        select(lmID,var,modifier) %>% 
                            arrange(lmID,var,modifier)
    
    l_var_compare <- setequal(df_var_1_exact,df_var_2_exact)
    res[['var']] <- if(l_var_compare==TRUE) TRUE else list(df_var_1_exact,df_var_2_exact)

#    message("variables compare: ", l_var_compare)

#3. compare lm_filter
    query_1 <- sprintf("SELECT * FROM %s WHERE lmID='%s' AND sessionID='%s'",tbl_lm_filter,lm_Id,session_Id_1)
    query_2 <- sprintf("SELECT * FROM %s WHERE lmID='%s' AND sessionID='%s'",tbl_lm_filter,lm_Id,session_Id_2)
    
    df_filt_1 <- dbGetQuery(db_conn, query_1)
    df_filt_2 <- dbGetQuery(db_conn, query_2)
    
    df_filt_1_exact <- df_filt_1 %>% select(lmID,filter_1,filter_2,filter_full)
    df_filt_2_exact <- df_filt_2 %>% select(lmID,filter_1,filter_2,filter_full)

    l_filt_compare <- setequal(df_filt_1_exact,df_filt_2_exact)
    res[['filt']] <- if(l_filt_compare==TRUE) TRUE else list(df_filt_1_exact,df_filt_2_exact)

#    message("filter compare: ", l_filt_compare)

#4. compare lm_interactions
    query_1 <- sprintf("SELECT * FROM %s WHERE lmID='%s' AND sessionID='%s'",tbl_lm_inter,lm_Id,session_Id_1)
    query_2 <- sprintf("SELECT * FROM %s WHERE lmID='%s' AND sessionID='%s'",tbl_lm_inter,lm_Id,session_Id_2)
    
    df_inter_1 <- dbGetQuery(db_conn, query_1)
    df_inter_2 <- dbGetQuery(db_conn, query_2)

    df_inter_1_exact <- df_inter_1 %>% select(lmID,var1,var2) %>% arrange(lmID,var1,var2)
    df_inter_2_exact <- df_inter_2 %>% select(lmID,var1,var2) %>% arrange(lmID,var1,var2)
    
    l_inter_compare <- setequal(df_inter_1_exact,df_inter_2_exact)
    res[['inter']] <- if(l_inter_compare==TRUE) TRUE else list(df_inter_1_exact,df_inter_2_exact)

#    message("interactions compare: ", l_inter_compare)

#5. compare lm_mutate
    query_1 <- sprintf("SELECT * FROM %s WHERE lmID='%s' AND sessionID='%s'",tbl_lm_mutate,lm_Id,session_Id_1)
    query_2 <- sprintf("SELECT * FROM %s WHERE lmID='%s' AND sessionID='%s'",tbl_lm_mutate,lm_Id,session_Id_2)
    
    df_mutate_1 <- dbGetQuery(db_conn, query_1)
    df_mutate_2 <- dbGetQuery(db_conn, query_2)
    
    df_mutate_1_exact <- df_mutate_1 %>% select(lmID,order,var,formula) %>% arrange(lmID,order,var,formula)
    df_mutate_2_exact <- df_mutate_2 %>% select(lmID,order,var,formula) %>% arrange(lmID,order,var,formula)
    
    l_mutate_compare <- setequal (df_mutate_1_exact, df_mutate_2_exact)
    res[['mutate']] <- if(l_mutate_compare==TRUE) TRUE else list(df_mutate_1_exact,df_mutate_2_exact)
    
#    message("mutate compare: ", l_mutate_compare)
    result[['lmID']]=lm_Id
    result[['site_1']]=site_Id_1
    result[['site_2']]=site_Id_2
    result[['session_ID_1']]=session_Id_1
    result[['session_ID_2']]=session_Id_2
    result[['ROI']]=ROI
    result[['result']]=l_main_lm_compare & l_var_compare & l_filt_compare & l_inter_compare & l_mutate_compare
    result[['details']]=res
    result
}
        
eregr_meta_int_select_latest_lm_results <- function(db_conn, study_Id,site_Id,lm_list,
                                                    tbl_res_keys="lm_results_keys", tbl_sess_lm="session_lm_analysis",
                                                    tbl_linear_model="linear_model") {
query <-sprintf("
SELECT DISTINCT LRK_MAX.lmID,LRK_MAX.metric,LRK_MAX.ROI, LRK_MAX.siteID, LRK_MAX.MAX_RSID, LRK_2.sessionID, LRK_2.result_sessionID,LDR.n_cont, LDR.n_pat, LDR.n_overall
FROM (SELECT LRK.lmID,LRK.metric,LRK.ROI, LM_SLA.siteID, MAX(LRK.result_sessionID) as MAX_RSID
		FROM (
			SELECT lmID,sessionID,SLA.studyID,siteID
			FROM (
				SELECT *
				FROM session_lm_analysis
				WHERE studyID='%s' AND siteID in ('%s') ) SLA, linear_model LM
				WHERE sessionID=SLA.session_analysis_ID AND lmID in ('%s')) LM_SLA, lm_results_keys LRK
		WHERE LM_SLA.lmID=LRK.lmID AND LM_SLA.sessionID=LRK.sessionID
		GROUP BY LRK.lmID,LRK.metric,LRK.ROI, LM_SLA.siteID ) LRK_MAX, lm_results_keys LRK_2, lm_demog_results LDR
WHERE LRK_MAX.MAX_RSID=LRK_2.result_sessionID 
AND LRK_MAX.lmID=LRK_2.lmID
AND LRK_MAX.metric=LRK_2.metric
AND LRK_MAX.ROI=LRK_2.ROI 
AND LRK_2.res_keyID=LDR.res_keyID",study_Id,paste(site_Id,collapse="','"),paste(lm_list,collapse="','"));
print(query)
dbGetQuery(db_conn,query)
}
eregr_meta_compare_model_across_sites <- function(db_conn, lm_Id, ROI_str, res_latest_lm) {
    res_ROI_LM<-res_latest_lm %>% 
            filter(ROI==ROI_str,lmID==lm_Id) %>% select (lmID,ROI,siteID,sessionID,MAX_RSID) %>% 
                distinct()
    res_ROI_LM$cj<-1
    r_cross<-inner_join(res_ROI_LM,res_ROI_LM,by='cj')%>% select(-cj) %>% filter(sessionID.x<sessionID.y)
    res<-pmap(list(lm_Id=r_cross$lmID.x,ROI=r_cross$ROI.x,site_Id_1=r_cross$siteID.x,site_Id_2=r_cross$siteID.y,
         session_Id_1=r_cross$sessionID.x,session_Id_2=r_cross$sessionID.y),eregr_meta_compare_mod_sess_pair,db_conn=db_conn)
}
eregr_meta_analysis_onemodel <- function(db_conn,study_Id,lm_Id,ROI_str,res_latest_lm, sess_id_and_time, compare_models=FALSE,
                                        tbl_session_meta="session_meta", tbl_meta_beta_res="meta_beta_results",
                                        tbl_meta_cohd_res="meta_cohd_results", tbl_sites_in_meta="sites_in_meta") {
    meta_onevertex <- function (df_vwise,metr_list,var_list,beta_name='beta',sterr_name='sterr') {         
        res_metrmap<-map(metr_list, function(cur_metr){
                res_varmap <- map(var_list,function(cur_var) {
                        cur_data <- df_vwise %>% filter (var==cur_var,metric==cur_metr)
                        cur_res=rma.uni(yi=cur_data[[beta_name]],sei=cur_data[[sterr_name]],slab=cur_data$siteID,method="REML",control=list(stepadj=0.5,maxiter=10000)) 
                        cur_res
                    })
            })
    }
    print(lm_Id)
    #compare model across sites - check that all sites have same version of model
   res_compare_str <- NA	 
   if(compare_models) {
	    res_compare_str <- eregr_meta_compare_model_across_sites(db_conn,lm_Id,ROI_str,res_latest_lm)
	    res_compare_lgl <- map_lgl(res_compare_str,"result")
#	    if(!all(res_compare_lgl)) stop(simpleError("not all models match")) #here we should return res_compare structure
    }
    #select compute meta-analysis for betas
    res_latest_filtered=res_latest_lm %>% 
                        filter(ROI==ROI_str,lmID==lm_Id) 
    res_session_ID <- res_latest_filtered %>%
                            select(result_sessionID) %>% 
                                distinct() %>% 
                                    unlist()
    query <- sprintf("
    SELECT res_keyID,lmID,siteID,ROI,metric,result_sessionID,vertex,var,beta,sterr 
    FROM lm_results LMR LEFT JOIN lm_results_keys USING (res_keyID) 
    LEFT JOIN session_lm_analysis ON lm_results_keys.sessionID=session_lm_analysis.session_analysis_ID
    WHERE result_sessionID in ('%s') AND lmID in ('%s') AND ROI='%s'
    ",paste(res_session_ID,collapse="','"),paste(lm_Id,collapse="','"),ROI_str)
    data <- dbGetQuery(db_conn,query)    #run rma.uni
    data_vwise_beta <- data %>% split (data$vertex)
    var_list <- unlist(select(data,var)%>%distinct())
    names(var_list) <- var_list
    metr_list <- unlist(select(data,metric)%>%distinct())
    names(metr_list) <- metr_list
    
    tic()
    data<-map(data_vwise_beta,safely(~meta_onevertex(.,metr_list,var_list)))
    data<-map(data,"result")
    data<-data[lapply(data,length)!=0]
    toc()

 #   return(data)
    #compute meta-analysis for cohen's d
        
    query <- sprintf("
    SELECT res_keyID,lmID,siteID,ROI,metric,result_sessionID,vertex,var,cohens_d,cohens_se
    FROM lm_cohend_results LMR LEFT JOIN lm_results_keys USING (res_keyID) 
    LEFT JOIN session_lm_analysis ON lm_results_keys.sessionID=session_lm_analysis.session_analysis_ID
    WHERE result_sessionID in ('%s') AND lmID in ('%s') AND ROI='%s'
    ",paste(res_session_ID,collapse="','"),paste(lm_Id,collapse="','"),ROI_str)
    data_cohd <- dbGetQuery(db_conn,query)    #run rma.uni
    data_vwise_cohd <- data_cohd %>% split (data_cohd$vertex)
    var_list_cohd <- unlist(select(data_cohd,var)%>%distinct())
    names(var_list_cohd) <- var_list_cohd
    metr_list_cohd <- unlist(select(data_cohd,metric)%>%distinct())
    names(metr_list_cohd) <- metr_list_cohd


    message("NOW TO COHEN's D!")
    tic()
    data_cohd<-map(data_vwise_cohd,safely(~ meta_onevertex(.,metr_list_cohd,var_list_cohd,'cohens_d','cohens_se')))
    data_cohd<-map(data_cohd,"result")
    data_cohd<-data_cohd[lapply(data_cohd,length)!=0]
   
    toc()
#    return(data_cohd)
    #write results to session_meta
    tbl_sess_meta <- eregr_int_get_tbl_header(db_conn,tbl_session_meta)
    tbl_sess_meta[1,'session_meta_ID']=sess_id_and_time[[2]]
    tbl_sess_meta[1,'studyID']=study_Id
    tbl_sess_meta[1,'lmID']=lm_Id
    tbl_sess_meta[1,'ROI']=ROI_str    
    print(tbl_sess_meta)
    l_sess_meta_write <- dbWriteTable(db_conn,name = tbl_session_meta,value = tbl_sess_meta,append=TRUE,row.names=FALSE)
    if(!l_sess_meta_write) stop("error writing meta-analysis session!")
    
    #write results to sites in meta-analysis
    query <- sprintf("
    SELECT * FROM %s 
    WHERE session_meta_ID='%s' AND studyID='%s' AND lmID='%s' AND ROI='%s'",
    tbl_session_meta,sess_id_and_time[[2]],study_Id,lm_Id,ROI_str)    
    tbl_sess_meta<-dbGetQuery(db_conn,query)
    
    tbl_sites <- tbl_sess_meta %>% 
                            inner_join(res_latest_filtered,by=c('lmID','ROI')) %>%
                                select(session_meta_ID,lmID,ROI,siteID,result_sessionID) %>%
                                    distinct()
                                
    l_sites_in_meta_write <- dbWriteTable(db_conn, name=tbl_sites_in_meta,value=tbl_sites,append=TRUE,row.names=FALSE)
    if(!l_sites_in_meta_write) stop ("error writing sites in metaanalysis table")
    l_sites_in_meta_write    

    #write results to meta_betas
    list_beta_res = list()
    for (i in seq_along(data)) {
        d_rb1<-pmap(list(d=data[[i]],n=names(data[[i]])), 
                                function(d,n){ 
                                res <- as.data.frame(do.call(rbind,d),stringsAsFactors=FALSE)
                                res[['var']] <- rownames(res)
                                res[['metric']] <- n
                                res[['vertex']] <- i
                                rownames(res) <- NULL
                                res})
        d_rb2 <- do.call(rbind,d_rb1)
        list_beta_res[[i]] <- d_rb2
    }
    df_beta_res <- as.data.frame(do.call(rbind,list_beta_res),stringsAsFactors=FALSE)
    if(nrow(df_beta_res)>0) {
	    df_beta_res[['meta_key_ID']] <- tbl_sess_meta[[1,'meta_key_ID']]
	    df_beta_res <- df_beta_res %>% 
        	                select(meta_key_ID,metric,vertex,var,b,se,zval,pval,ci.lb,ci.ub,tau2,se.tau2,I2,H2) %>%
                	            rename(ci_lb = ci.lb,ci_ub=ci.ub,se_tau2=se.tau2)
	    df_beta_res[,5:14]<-sapply(df_beta_res[,5:14],as.numeric)
	    rownames(df_beta_res) <- NULL
#    return(df_beta_res)
	    l_write_beta_res <- dbWriteTable(db_conn,name=tbl_meta_beta_res,value = df_beta_res,append=TRUE,row.names=FALSE)
    }
    else
	l_write_beta_res<-FALSE
    if(!l_write_beta_res) stop("error writing beta - metaanalysis results")
    l_write_beta_res
    #write results to meta_cohend
    list_cohd_res = list()
    for (i in seq_along(data_cohd)) {
        d_rb1<-pmap(list(d=data_cohd[[i]],n=names(data_cohd[[i]])), 
                                function(d,n){ 
                                res <- as.data.frame(do.call(rbind,d),stringsAsFactors=FALSE)
                                res[['var']] <- rownames(res)
                                res[['metric']] <- n
                                res[['vertex']] <- i
                                rownames(res) <- NULL
                                res})
        d_rb2 <- do.call(rbind,d_rb1)
        list_beta_res[[i]] <- d_rb2
    }
    df_cohd_res <- as.data.frame(do.call(rbind,list_beta_res),stringsAsFactors=FALSE)
    if(nrow(df_cohd_res)>0) {
	    df_cohd_res[['meta_key_ID']] <- tbl_sess_meta[[1,'meta_key_ID']]
	    df_cohd_res <- df_cohd_res %>% 
        	                select(meta_key_ID,metric,vertex,var,b,se,zval,pval,ci.lb,ci.ub,tau2,se.tau2,I2,H2) %>%
                	            rename(cohd=b,ci_lb = ci.lb,ci_ub=ci.ub,se_tau2=se.tau2)
	    df_cohd_res[,5:14]<-sapply(df_cohd_res[,5:14],as.numeric)
	    rownames(df_cohd_res) <- NULL
	#    return(df_beta_res)
	    l_write_cohd_res <- dbWriteTable(db_conn,name=tbl_meta_cohd_res,value = df_cohd_res,append=TRUE,row.names=FALSE)        
    }
    else 
	l_write_cohd_res <- FALSE
    list(l_beta=l_write_beta_res,l_cohd=l_write_cohd_res, res_compare = res_compare_str)
    
#    res_compare_lgl
}

#forest plots
demog_summary <- function (db_conn, study_Id, lm_Id, ROI,metric,site_list, use_meta=FALSE) {
    query <- sprintf("SELECT LRK.res_keyID,vertex,var,SLA.studyID,SLA.siteID,LM.lmID,LDR.n_cont, LDR.n_pat, LDR.n_overall, cohens_d,cohens_se,cohens_low_ci,cohens_high_ci,cohens_pval 
FROM lm_cohend_results LCR, lm_demog_results LDR, lm_results_keys LRK, linear_model LM, session_lm_analysis SLA, (SELECT studyID,siteID,lmID,ROI,metric,MAX(result_sessionID) as max_res_sess_ID FROM lm_results_keys LRK, session_lm_analysis SLA
WHERE LRK.sessionID=SLA.session_analysis_ID
AND SLA.studyID='%s' AND LRK.lmID='%s'
AND ROI='%s'
AND metric='%s'
AND SLA.siteID IN ('%s')
GROUP BY studyID,siteID,lmID,ROI,metric) RES_MAX
WHERE LCR.res_keyID=LRK.res_keyID AND LDR.res_keyID=LRK.res_keyID AND LDR.res_keyID=LCR.res_keyID AND LRK.lmID=LM.lmID and LRK.sessionID=SLA.session_analysis_ID and 
LRK.sessionID=LM.sessionID
AND SLA.studyID='%s' and LM.lmID='%s'
AND LRK.result_sessionID = max_res_sess_ID 
AND SLA.siteID = RES_MAX.siteID and LRK.lmID=RES_MAX.lmID and LRK.ROI=RES_MAX.ROI and LRK.metric=RES_MAX.metric
AND SLA.studyID=RES_MAX.studyID",study_Id,lm_Id,ROI,metric,paste(site_list,collapse="','"),study_Id,lm_Id)
    res <- dbGetQuery(db_conn,query) %>% arrange(siteID)



}


forest_plot_cohd <- function(db_conn, study_Id, lm_Id, ROI,metric,site_list, use_meta=FALSE,rename_site_map=NA) {
    message(lm_Id)
    query <- sprintf("SELECT LRK.res_keyID,vertex,var,SLA.studyID,SLA.siteID,LM.lmID,cohens_d,cohens_se,cohens_low_ci,cohens_high_ci,cohens_pval 
FROM lm_cohend_results LCR, lm_results_keys LRK, linear_model LM, session_lm_analysis SLA, (SELECT studyID,siteID,lmID,ROI,metric,MAX(result_sessionID) as max_res_sess_ID FROM lm_results_keys LRK, session_lm_analysis SLA
WHERE LRK.sessionID=SLA.session_analysis_ID
AND SLA.studyID='%s' AND LRK.lmID='%s'
AND ROI='%s'
AND metric='%s'
AND SLA.siteID IN ('%s')
GROUP BY studyID,siteID,lmID,ROI,metric) RES_MAX
WHERE LCR.res_keyID=LRK.res_keyID AND LRK.lmID=LM.lmID and LRK.sessionID=SLA.session_analysis_ID and 
LRK.sessionID=LM.sessionID
AND SLA.studyID='%s' and LM.lmID='%s'
AND LRK.result_sessionID = max_res_sess_ID 
AND SLA.siteID = RES_MAX.siteID and LRK.lmID=RES_MAX.lmID and LRK.ROI=RES_MAX.ROI and LRK.metric=RES_MAX.metric
AND SLA.studyID=RES_MAX.studyID",study_Id,lm_Id,ROI,metric,paste(site_list,collapse="','"),study_Id,lm_Id)
    res <- dbGetQuery(db_conn,query) %>% arrange(siteID)
print(query)
print(res)
print(rename_site_map)
    if(!is.na(rename_site_map)) {
	rename_map=read.csv(rename_site_map)
	print(rename_map)
	renamer <-function(site_id,match_id) {
		if (is.na(match_id)){
			return (site_id)
		}
		else 
			print(match_id)
print(rename_map$new[match_id])
			return (as.character(rename_map$new[match_id]))

	}
	res$siteID=map2_chr(res$siteID,match(res$siteID,rename_map$old),~renamer(.x,.y))
	print(res$siteID)
    }
    if(use_meta) {
        query_meta <- sprintf("SELECT * FROM meta_cohd_results MCR, session_meta SM, (SELECT MAX(session_meta_ID) AS max_meta_ID FROM session_meta SM 
WHERE studyID='%s' AND lmID='%s' AND ROI='%s'
) SM_MAX
WHERE MCR.meta_key_ID=SM.meta_key_ID AND 
SM.session_meta_ID=max_meta_ID
AND studyID='%s' AND lmID='%s' AND ROI='%s' AND metric='%s'",study_Id,lm_Id,ROI,study_Id,lm_Id,ROI,metric)
        meta_res <- dbGetQuery(db_conn,query_meta)
        return (forestplot(title = sprintf( "Cohen's D summary\n ROI: %s\n metric: %s\n model: %s",ROI,metric,lm_Id),labeltext = c(res$siteID,"SUMMARY"),
           mean = c(res$cohens_d,meta_res$cohd),
           lower  = c(res$cohens_low_ci,meta_res$ci_lb), 
           upper = c(res$cohens_high_ci,meta_res$ci_ub),
           is.summary = c(rep(FALSE,length(res$cohens_d)),rep(TRUE,length(meta_res$cohd))),
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue")))
    }
        
    forestplot(title = sprintf("Cohen's D summary\n ROI: %s\n metric: %s\n model: %s",ROI,metric,lm_Id),
           labeltext = c(res$siteID),
           mean = c(res$cohens_d),
           lower  = c(res$cohens_low_ci), 
           upper = c(res$cohens_high_ci),
           is.summary = c(rep(FALSE,length(res$cohens_d))),
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"))
}



#feature set functionality
eregr_register_feature_sets <- function (db_conn, study_Id,site_Id, 
                                         erase_existing_fs = TRUE, 
                                         tbl_fs_name="feature_sets",tbl_study_metr='study_metrics', 
                                         tbl_sites_in_study ='sites_in_study', tbl_session_fs='session_fs',
                                        tbl_fs_newregr='feature_sets_newregr') {
    study_row <- eregr_get_study_info(db_conn,study_Id)
#    metric_list <- dbReadTable(db_conn,name=tbl_study_metr) %>%
#                            filter (studyID==study_Id) %>%
#                                    select(metr_name) %>%
#                                            unlist()

# read feature set gdoc
    fs_path <- study_row[['study_fs_path']]
    fs_data <- fs_path %>%
                    gs_url() %>%
                            gs_read()
# get study_site_ID
    query <- sprintf("SELECT study_site_ID FROM %s
         WHERE studyID='%s' AND siteID='%s'",tbl_sites_in_study,study_Id,site_Id)
    study_site_ID <- dbGetQuery(db_conn,query)
    if (nrow(study_site_ID)<1) 
        stop ("couldn't find study_site_ID for study ", study_Id, " and site ", site_Id)
    else if (nrow(study_site_ID)>1) 
         stop ("Found more than 1 study_site_ID for study ", study_Id, " and site ", site_Id ,". CHECK DATABASE CONSISTENCY!")
    study_site_ID <- study_site_ID[[1]]

# generate session_fs_ID:
    time_and_id <- eregr_int_get_unique_time_id(db_conn)
    session_fs_ID <- time_and_id[[2]]
    timestamp <- time_and_id[[1]]

    apply(fs_data,1, function(fs_line) {
        if (is.na(fs_line[['FeatureSet']])) 
            return (NA)
        fsID <- fs_line[['FeatureSet']]
        
        #preparing and writing session into session_fs table
        session_fs_tbl <- eregr_int_get_tbl_header(db_conn,tbl_session_fs)
        session_fs_tbl[1,'session_fs_ID'] <- session_fs_ID
        session_fs_tbl[1,'fsID'] <- fsID
        session_fs_tbl[1,'timestamp'] <- timestamp
        session_fs_tbl[1,'study_site_ID'] <- study_site_ID
                    
        l_session_fs <- dbWriteTable(db_conn, name=tbl_session_fs, value = session_fs_tbl, append = TRUE, row.names = FALSE)
#        l_session_fs
        
        # writing combinations of ROI and metrics for feature set
        fs_line[['Contents']] <- str_replace_all(fs_line[['Contents']],pattern='\\s',replacement="")
	fs_line[['Metrics']] <- str_replace_all(fs_line[['Metrics']],pattern='\\s',replacement="")
	metric_list <- strsplit(fs_line[['Metrics']],split=',')[[1]]
	print(metric_list)

    	if (length(metric_list) <1 ) stop("length of metrics_list for study ", study_Id, " less than 1")

	ROI_list <- strsplit(fs_line[['Contents']],split=',')[[1]]
        #check for ROI length > 0
        if(length(ROI_list)<1 | sum(is.na(ROI_list))>0) 
            stop("could not extract ROI list for feature set: ", fsID)

        l_reg_ROI <- eregr_register_ROI_metr_for_FS(db_conn,fsID,session_fs_ID,metric_list,ROI_list)
            
        #preparing and writing new regressors into feature_sets_newregr table
        
        newregr <- str_replace_all(fs_line[['NewRegressors']],pattern='\\s',replacement="")
        if(is.na(newregr) | newregr=="") return(NA)

        split_regr <- str_split(newregr,';')[[1]]
        command_to_exec <- str_replace_all(split_regr,"^.*=","")
        command_to_exec <- str_replace_all(command_to_exec,"(?:\\{ROI:|\\{cov:)([\\w]+)\\}",str_c("as.numeric(","\\1",")"))
        new_regr_name <- str_replace_all(split_regr,"=.*$","") 
        
            
        fs_newregr_tbl <- eregr_int_get_tbl_header(db_conn, tbl_fs_newregr)
        fs_newregr_tbl[1:length(new_regr_name),'fsID'] <- NA
        fs_newregr_tbl[['var']] <- new_regr_name
        fs_newregr_tbl[['formula']] <- command_to_exec
        fs_newregr_tbl[['fsID']] <- fsID
        fs_newregr_tbl[['session_fs_ID']] <- session_fs_ID
        l_fs_newregr <- dbWriteTable(db_conn,name=tbl_fs_newregr,value = fs_newregr_tbl, append = TRUE, row.names = FALSE)

	fs_newregr_tbl        
#        list(l_session_fs, l_reg_ROI, l_fs_newregr)
    })
}


eregr_int_compute_new_covariates <- function (db_conn, newregr_list, study_Id, site_Id, tbl_fs_covars = "feature_sets_covariates") {
    df_covars <- eregr_get_covariates(db_conn,study_Id,site_Id)
    df_covars <- df_covars %>% spread(cov_name,cov_value)

    df_metr <- eregr_get_all_metrics_data(db_conn,study_Id,site_Id)
    uniq_metr <- unique(df_metr[['metric']])
    uniq_vertex <- unique (df_metr[['vertex']])
    fs_list <- list()
    if(length(uniq_vertex)>1) stop("as of now, feature sets work only with ROI analysis, shape analysis is not supported")
    for (elem in uniq_metr) {
        df_cur_metr <- df_metr %>% 
                        filter (metric==elem) %>% 
                            select(subjID,ROI,value) %>%
                                spread(ROI,value)
        df_metr_covar <- inner_join(df_cur_metr,df_covars,by='subjID')        

        fs_list[[elem]] <- map(newregr_list, function(df,df_newregr_tbl) {
                if(class(df)!="data.frame")
                   return (NA)                
                for (i in 1:nrow(df)) {
                    cmd <- str_c("mutate(df_newregr_tbl,",df[[i,'var']],"=",df[[i,'formula']],")")
                    df_newregr_tbl <- eval(parse(text=cmd))
                }
                #extract only new regressors values
                df_varstr <-paste(df[['var']],collapse=',')
                cmd=paste("select(df_newregr_tbl,subjID,",df_varstr,")",sep="")                
                df_newregr_tbl <- eval(parse(text=cmd))
                df_newregr_tbl <- df_newregr_tbl %>% gather("cov_name","cov_value",2:ncol(df_newregr_tbl))
                fsID <- df[[1,'fsID']]
                session_fs_ID <- df [[1,'session_fs_ID']]
                    
                df_newregr_tbl[['fsID']] <- fsID
                df_newregr_tbl[['session_fs_ID']] <- session_fs_ID
                df_newregr_tbl[['metric']] <- elem
                df_newregr_tbl[['vertex']] <- 1 #as of now we support only 1-vertex ROI analysis
                df_newregr_tbl
                },df_newregr_tbl = df_metr_covar)
        
    }
    # write to database
    l_fs_cov = TRUE
    dbBegin(db_conn)
    for (elem in uniq_metr) {
        for (i in seq_along(fs_list[[elem]])) {
            if (!( class(fs_list[[elem]][[i]])=="data.frame"))
                next
            rec <- fs_list[[elem]][[i]]
            rec[['row_ID']]<-NA
            rec <- rec %>% select(session_fs_ID,fsID,subjID,cov_name,cov_value,metric,vertex)
            
                #row_ID set to NA
            l_fs_cov <- l_fs_cov &  dbWriteTable(db_conn, name=tbl_fs_covars, value = rec, append=TRUE, row.names=FALSE)
           
        }
    }
    if (l_fs_cov) {
		message("commiting transactions")
    		dbCommit(db_conn)
	}
    else {
	message("transcations rollback")
	dbRollback(db_conn)
    }
	
    fs_list                   
}



eregr_register_ROI_metr_for_FS <- function (db_conn, fs_Id, session_fs_ID, metric_list,ROI_list,tbl_fs_name="feature_sets") {
        tbl_fs <- eregr_int_get_tbl_header(db_conn,tbl_fs_name)
        c = 0
    	for (j_metr in seq_along(metric_list))
                for(k_ROI in seq_along(ROI_list)) {
                    c=c+1
                    tbl_fs[c,'metric']=metric_list[[j_metr]]
                    tbl_fs[c,'ROI']=ROI_list[[k_ROI]]
                    tbl_fs[c,'session_fs_ID']=session_fs_ID
                    tbl_fs[c,'fsID']=fs_Id
                }
        l_write_tbl_fs <- dbWriteTable(db_conn,name = tbl_fs_name, value = tbl_fs, append=TRUE, row.names=FALSE)     
        l_write_tbl_fs
}

eregr_get_fs_for_lm_ROI <- function (db_conn, fs_Id,session_fs_Id,ROI,tbl_fs_name="feature_sets") {
    query <- sprintf("SELECT * 
                      FROM %s FS
                      WHERE fsID='%s' AND session_fs_ID = '%s' AND ROI = '%s' ",tbl_fs_name,fs_Id, session_fs_Id, ROI)
    dbGetQuery(db_conn,query)        
}
    
