library(shiny)
library(datasets)
library(RMySQL)
library(dplyr)
library(tidyr)

# Shiny app to export regression coefficients from meta-analysis stage
# Dmitry Isaev
# Imaging Genetics Center, Institute of Neuroimaging and Informatics (INI), 
# University of Southern California
# 2018

shinyServer(function(input, output) {
  
  #connect to database
  prepare_meta_beta_summ <- function() {
    # !!! replace this path to path to your 'lib_eregr_core.R'
    suppressPackageStartupMessages(source(file='/path_to_lib_eregr_core/lib_eregr_core.R'))
    # !!! replace this path to path to your 'database_connect.R'
    suppressPackageStartupMessages(source(file='/path_to_database_connect/database_connect_mdd.R'))
    
    db_conn <- database_connect()
    on.exit(dbDisconnect(db_conn))
    
    query <- sprintf("
                  SELECT * FROM session_meta SM,
			                        meta_beta_results MBR,
                     ( SELECT MAX_SM.max_smID,LRK.lmID,LRK.ROI,LRK.metric,SUM(n_cont),SUM(n_pat),SUM(n_overall) FROM sites_in_meta SIS,
                     session_meta SM,
                     ( SELECT MAX(session_meta_ID) as max_smID FROM session_meta
                     WHERE lmID='%s' AND studyID='%s') MAX_SM,
                     lm_results_keys LRK,
                     lm_demog_results LDR
                     WHERE SIS.session_meta_ID=MAX_SM.max_smID
                     AND SIS.lmID='%s'
                     AND LRK.result_sessionID=SIS.result_sessionID
                     AND LRK.lmID=SIS.lmID
                     AND LRK.ROI=SIS.ROI
                     AND LRK.metric='%s'
                     #AND LRK.ROI='ACR'
                     AND LDR.res_keyID=LRK.res_keyID
                     AND SM.session_meta_ID=MAX_SM.max_smID
                     AND SM.lmID=SIS.lmID
                     AND SM.ROI=SIS.ROI
                     GROUP BY MAX_SM.max_smID,LRK.lmID, LRK.ROI, LRK.metric
                     ) N_DEMOG
                     WHERE SM.session_meta_ID=N_DEMOG.max_smID
                     AND SM.ROI=N_DEMOG.ROI
                     AND SM.lmID=N_DEMOG.lmID
                     AND MBR.metric=N_DEMOG.metric
                     AND MBR.meta_key_ID=SM.meta_key_ID
                     ",input$lm,input$studyID,input$lm,input$metric)
    res <- dbGetQuery(db_conn,query)
    res
}  
  output$meta_beta_summary <- renderTable({
    prepare_meta_beta_summ()
    
  },digits = 8)
})