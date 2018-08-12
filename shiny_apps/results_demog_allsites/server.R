library(shiny)
library(datasets)
library(RMySQL)
library(dplyr)
library(tidyr)

# Shiny app - Quality control of amount of patients/controls in each ROI in the same model
# Dmitry Isaev
# Imaging Genetics Center, Institute of Neuroimaging and Informatics (INI), 
# University of Southern California
# 2018

shinyServer(function(input, output) {
  
  prepare_demog_summ <- function() {
    # !!! replace this path to path to your 'lib_eregr_core.R'
    suppressPackageStartupMessages(source(file='/path_to_lib_eregr_core/lib_eregr_core.R'))
    # !!! replace this path to path to your 'database_connect.R'
    suppressPackageStartupMessages(source(file='/path_to_database_connect/database_connect_mdd.R'))
    
    db_conn <- database_connect()
    on.exit(dbDisconnect(db_conn))
    
    query <- sprintf("
SELECT * FROM (SELECT lmID,siteID,ROI,MAX(result_sessionID) max_rsID from lm_results_keys LRK2,session_lm_analysis SLA2
	WHERE LRK2.lmID='%s'
                     AND  LRK2.sessionID=SLA2.session_analysis_ID
                     AND SLA2.studyID='%s'
                     GROUP BY lmID,siteID,ROI) MAX_RS,
                     lm_results_keys LRK,
                     lm_demog_results LDR
                     WHERE LRK.result_sessionID=MAX_RS.max_rsID
                     AND LRK.lmID=MAX_RS.lmID
                     AND LRK.ROI=MAX_RS.ROI
                     AND LDR.res_keyID=LRK.res_keyID
                     AND siteID='%s'",input$lm,input$studyID,input$site)
    
    res <- dbGetQuery(db_conn,query)
    res
}  
  output$demog_summary <- renderTable({
    prepare_demog_summ()
    
  },digits = 8)
})