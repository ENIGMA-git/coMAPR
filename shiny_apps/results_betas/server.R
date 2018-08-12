library(shiny)
library(datasets)
library(RMySQL)
library(dplyr)
library(tidyr)

# Shiny app to output regression coefficients with standard errors
# Dmitry Isaev
# Imaging Genetics Center, Institute of Neuroimaging and Informatics (INI), 
# University of Southern California
# 2018

shinyServer(function(input, output) {
  
  prepare_betas_summ <- function() {
    # !!! replace this path to path to your 'lib_eregr_core.R'
    suppressPackageStartupMessages(source(file='/path_to_lib_eregr_core/lib_eregr_core.R'))
    # !!! replace this path to path to your 'database_connect.R'
    suppressPackageStartupMessages(source(file='/path_to_database_connect/database_connect.R'))
    
    db_conn <- database_connect()
    on.exit(dbDisconnect(db_conn))
    
    query <- sprintf("SELECT SLA.studyID,SLA.siteID,LRK.res_keyID,LRK.result_sessionID,var,beta,sterr,p_beta FROM lm_results_keys LRK,
(SELECT MAX(result_sessionID) as max_res_sesID FROM lm_results_keys,session_lm_analysis
WHERE lmID='%s' AND studyID='%s' AND siteID='%s' AND sessionID=session_analysis_ID AND ROI='%s') MAX_LRK,lm_results LR,session_lm_analysis SLA
WHERE LRK.result_sessionID=MAX_LRK.max_res_sesID AND lmID='%s'
AND metric='%s' AND ROI='%s'
AND LRK.res_keyID=LR.res_keyID AND LRK.sessionID=SLA.session_analysis_ID",input$lm,input$studyID,input$site,input$roi,input$lm, input$metric,input$roi)
    
    res <- dbGetQuery(db_conn,query)
    res
}  
  output$betas_summary <- renderTable({
    prepare_betas_summ()
    
  },digits = 8)
})
