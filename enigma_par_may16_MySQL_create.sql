-- Created by Vertabelo (http://vertabelo.com)
-- Last modification date: 2017-05-25 23:33:34.768

-- tables
-- Table: covariates_general
CREATE TABLE covariates_general (
    subjID character(100) NOT NULL,
    cov_name character(100) NOT NULL,
    cov_value character(255) NULL,
    session_covar_ID character(200) NOT NULL,
    CONSTRAINT covariates_general_pk PRIMARY KEY (subjID,cov_name,session_covar_ID)
);

CREATE INDEX covariates_general_idx_1 ON covariates_general (subjID,cov_name);

-- Table: feature_sets
CREATE TABLE feature_sets (
    fsID character(100) NOT NULL,
    studyID character(100) NOT NULL,
    lmID character(100) NOT NULL,
    metric character(100) NOT NULL,
    ROI character(100) NOT NULL,
    CONSTRAINT feature_sets_pk PRIMARY KEY (fsID,studyID,lmID,metric,ROI)
);

-- Table: linear_model
CREATE TABLE linear_model (
    lmID character(100) NOT NULL,
    sessionID character(200) NOT NULL,
    lm_gDoc_path varchar(1000) NOT NULL,
    lm_name varchar(1000) NOT NULL,
    lm_text varchar(1000) NOT NULL,
    main_factor character(100) NULL,
    new_regressors varchar(1000) NOT NULL,
    cont_value integer NULL DEFAULT 0,
    pat_value integer NULL DEFAULT 1,
    cont_min integer NULL DEFAULT 0,
    pat_min integer NULL DEFAULT 0,
    comments varchar(1000) NULL,
    CONSTRAINT linear_model_pk PRIMARY KEY (lmID,sessionID)
);

CREATE INDEX linear_model_idx_1 ON linear_model (sessionID);

CREATE INDEX linear_model_idx_2 ON linear_model (lmID,sessionID);

-- Table: lm_cohend_results
CREATE TABLE lm_cohend_results (
    res_keyID bigint NOT NULL,
    vertex integer NOT NULL,
    var character(100) NOT NULL,
    cohens_d double NOT NULL,
    cohens_se double NOT NULL,
    cohens_low_ci double NOT NULL,
    cohens_high_ci double NOT NULL,
    cohens_pval double NOT NULL,
    CONSTRAINT lm_cohend_results_pk PRIMARY KEY (res_keyID,vertex,var)
);

-- Table: lm_corr_results
CREATE TABLE lm_corr_results (
    res_keyID bigint NOT NULL,
    vertex character(100) NOT NULL,
    var character(100) NOT NULL,
    corr double NOT NULL,
    corr_pval double NOT NULL,
    corr_se double NOT NULL,
    CONSTRAINT lm_corr_results_pk PRIMARY KEY (res_keyID,vertex,var)
);

-- Table: lm_demog_results
CREATE TABLE lm_demog_results (
    res_keyID bigint NOT NULL,
    n_overall integer NOT NULL,
    n_cont integer NOT NULL,
    n_pat integer NOT NULL,
    CONSTRAINT lm_demog_results_pk PRIMARY KEY (res_keyID)
);

-- Table: lm_filter
CREATE TABLE lm_filter (
    lmID character(100) NOT NULL,
    sessionID character(200) NOT NULL,
    filter_1 varchar(1000) NULL,
    filter_2 varchar(1000) NULL,
    filter_full varchar(1000) NULL,
    CONSTRAINT lm_filter_pk PRIMARY KEY (lmID,sessionID)
);

CREATE INDEX lm_filter_idx_1 ON lm_filter (lmID,sessionID);

-- Table: lm_interactions
CREATE TABLE lm_interactions (
    var1 character(100) NOT NULL,
    var2 character(100) NOT NULL,
    lmID character(100) NOT NULL,
    sessionID character(100) NOT NULL,
    CONSTRAINT lm_interactions_pk PRIMARY KEY (var1,var2,lmID,sessionID)
);

CREATE INDEX lm_interactions_idx_1 ON lm_interactions (lmID,sessionID);

-- Table: lm_mutate
CREATE TABLE lm_mutate (
    formula varchar(1000) NOT NULL,
    `order` integer NOT NULL,
    var character(100) NOT NULL,
    lmID character(100) NOT NULL,
    sessionID character(200) NOT NULL,
    CONSTRAINT lm_mutate_pk PRIMARY KEY (var,lmID,sessionID)
);

CREATE INDEX lm_mutate_idx_1 ON lm_mutate (lmID,sessionID);

-- Table: lm_results
CREATE TABLE lm_results (
    res_keyID bigint NOT NULL,
    vertex integer NOT NULL,
    var character(100) NOT NULL,
    beta double NOT NULL,
    sterr double NOT NULL,
    p_beta double NOT NULL,
    CONSTRAINT lm_results_pk PRIMARY KEY (res_keyID,vertex,var)
);

-- Table: lm_results_keys
CREATE TABLE lm_results_keys (
    res_keyID bigint NOT NULL AUTO_INCREMENT,
    lmID character(100) NOT NULL,
    sessionID character(200) NOT NULL,
    metric char(100) NOT NULL,
    ROI character(100) NOT NULL,
    result_sessionID character(200) NOT NULL,
    CONSTRAINT lm_results_keys_pk PRIMARY KEY (res_keyID)
);

-- Table: lm_summary_cont
CREATE TABLE lm_summary_cont (
    res_keyID bigint NOT NULL,
    var character(100) NOT NULL,
    mean double NOT NULL,
    median double NOT NULL,
    Q1 double NOT NULL,
    Q3 double NOT NULL,
    min double NOT NULL,
    max double NOT NULL,
    std double NOT NULL,
    CONSTRAINT lm_summary_cont_pk PRIMARY KEY (res_keyID,var)
);

-- Table: lm_summary_fact
CREATE TABLE lm_summary_fact (
    res_keyID bigint NOT NULL,
    var character(100) NOT NULL,
    value int NOT NULL,
    amount int NOT NULL,
    CONSTRAINT lm_summary_fact_pk PRIMARY KEY (res_keyID,var,value)
);

-- Table: lm_variables
CREATE TABLE lm_variables (
    var character(100) NOT NULL,
    lmID character(100) NOT NULL,
    sessionID character(200) NOT NULL,
    modifier character(100) NULL,
    CONSTRAINT lm_variables_pk PRIMARY KEY (var,lmID,sessionID)
);

CREATE INDEX lm_variables_idx_1 ON lm_variables (lmID,sessionID);

-- Table: meta_beta_results
CREATE TABLE meta_beta_results (
    meta_key_ID int NOT NULL,
    metric char(100) NOT NULL,
    vertex int NOT NULL,
    var character(100) NOT NULL,
    b double NOT NULL,
    se double NOT NULL,
    zval double NOT NULL,
    pval double NOT NULL,
    ci_lb double NOT NULL,
    ci_ub double NOT NULL,
    tau2 double NOT NULL,
    se_tau2 double NOT NULL,
    I2 double NOT NULL,
    H2 double NOT NULL,
    CONSTRAINT meta_beta_results_pk PRIMARY KEY (meta_key_ID,metric,vertex,var)
);

-- Table: meta_cohd_results
CREATE TABLE meta_cohd_results (
    meta_key_ID int NOT NULL,
    metric char(100) NOT NULL,
    vertex int NOT NULL,
    var character(100) NOT NULL,
    cohd double NOT NULL,
    se double NOT NULL,
    zval double NOT NULL,
    pval double NOT NULL,
    ci_lb double NOT NULL,
    ci_ub double NOT NULL,
    tau2 double NOT NULL,
    se_tau2 double NOT NULL,
    I2 double NOT NULL,
    H2 double NOT NULL,
    CONSTRAINT meta_cohd_results_pk PRIMARY KEY (meta_key_ID,metric,vertex,var)
);

-- Table: metrics_data
CREATE TABLE metrics_data (
    session_metr_ID character(200) NOT NULL,
    subjID character(100) NOT NULL,
    metric character(100) NOT NULL,
    ROI character(100) NOT NULL,
    vertex integer NOT NULL,
    value double NULL,
    row_ID integer NOT NULL AUTO_INCREMENT,
    CONSTRAINT metrics_data_pk PRIMARY KEY (session_metr_ID,row_ID)
);

CREATE INDEX metrics_data_idx_1 ON metrics_data (subjID,metric,ROI);

CREATE INDEX metrics_data_idx_2 ON metrics_data (row_ID);

-- Table: metrics_file
CREATE TABLE metrics_file (
    studyID character(100) NOT NULL,
    siteID character(100) NOT NULL,
    subjID character(100) NOT NULL,
    metric character(100) NOT NULL,
    file_path varchar(1000) NOT NULL,
    file_lastedit datetime NOT NULL,
    ROI character(100) NOT NULL,
    session_metr_ID character(200) NOT NULL,
    CONSTRAINT metrics_file_pk PRIMARY KEY (studyID,siteID,subjID,metric,ROI)
);

CREATE INDEX metrics_file_idx_1 ON metrics_file (studyID,siteID,subjID,ROI);

-- Table: session_covariates
CREATE TABLE session_covariates (
    session_covar_ID character(200) NOT NULL,
    timestamp datetime NOT NULL,
    file_cov_path varchar(1000) NOT NULL,
    file_lastedit datetime NOT NULL,
    study_site_ID int NOT NULL,
    CONSTRAINT session_covariates_pk PRIMARY KEY (session_covar_ID)
);

CREATE INDEX session_covariates_idx_1 ON session_covariates (session_covar_ID);

-- Table: session_lm_analysis
CREATE TABLE session_lm_analysis (
    session_analysis_ID character(200) NOT NULL,
    studyID character(100) NOT NULL,
    siteID character(100) NOT NULL,
    timestamp datetime NOT NULL,
    CONSTRAINT session_lm_analysis_pk PRIMARY KEY (session_analysis_ID)
);

CREATE INDEX session_lm_analysis_idx_1 ON session_lm_analysis (session_analysis_ID);

-- Table: session_loadMetrics
CREATE TABLE session_loadMetrics (
    session_metr_ID character(200) NOT NULL,
    timestamp datetime NOT NULL,
    study_site_ID int NOT NULL,
    CONSTRAINT session_loadMetrics_pk PRIMARY KEY (session_metr_ID)
);

-- Table: session_meta
CREATE TABLE session_meta (
    meta_key_ID int NOT NULL AUTO_INCREMENT,
    session_meta_ID character(200) NOT NULL,
    studyID character(100) NOT NULL,
    lmID character(100) NOT NULL,
    ROI character(100) NOT NULL,
    UNIQUE INDEX session_meta_ak_1 (session_meta_ID,lmID,ROI),
    CONSTRAINT session_meta_pk PRIMARY KEY (meta_key_ID)
);

-- Table: sites_in_meta
CREATE TABLE sites_in_meta (
    session_meta_ID character(200) NOT NULL,
    lmID character(100) NOT NULL,
    ROI character(100) NOT NULL,
    siteID character(100) NOT NULL,
    result_sessionID character(200) NOT NULL,
    CONSTRAINT sites_in_meta_pk PRIMARY KEY (session_meta_ID,lmID,ROI,siteID)
);

-- Table: sites_in_study
CREATE TABLE sites_in_study (
    siteID character(100) NOT NULL,
    studyID character(100) NOT NULL,
    study_site_ID int NULL AUTO_INCREMENT,
    UNIQUE INDEX study_siteID (study_site_ID),
    CONSTRAINT sites_in_study_pk PRIMARY KEY (siteID,studyID)
);

-- Table: studies
CREATE TABLE studies (
    studyID character(100) NOT NULL,
    study_name varchar(100) NULL,
    study_gDoc_path varchar(1000) NOT NULL,
    study_analysis_path varchar(1000) NOT NULL,
    study_demographics_path varchar(1000) NOT NULL,
    study_data_format character(100) NOT NULL,
    CONSTRAINT studies_pk PRIMARY KEY (studyID)
);

CREATE INDEX studies_idx_1 ON studies (studyID);

-- Table: study_metrics
CREATE TABLE study_metrics (
    studyID character(100) NOT NULL,
    metr_name character(100) NOT NULL,
    CONSTRAINT study_metrics_pk PRIMARY KEY (studyID,metr_name)
);

CREATE INDEX study_metrics_idx_1 ON study_metrics (studyID,metr_name);

-- Table: subj_excluded
CREATE TABLE subj_excluded (
    siteID character(100) NOT NULL,
    studyID character(100) NOT NULL,
    subjID character(100) NOT NULL,
    ROI character(100) NOT NULL,
    QC integer NULL,
    session_metr_ID character(200) NOT NULL,
    CONSTRAINT subj_excluded_pk PRIMARY KEY (siteID,studyID,subjID,ROI)
);

-- foreign keys
-- Reference: Table_14_studies (table: study_metrics)
ALTER TABLE study_metrics ADD CONSTRAINT Table_14_studies FOREIGN KEY Table_14_studies (studyID)
    REFERENCES studies (studyID);

-- Reference: covariates_general_session_covariates (table: covariates_general)
ALTER TABLE covariates_general ADD CONSTRAINT covariates_general_session_covariates FOREIGN KEY covariates_general_session_covariates (session_covar_ID)
    REFERENCES session_covariates (session_covar_ID);

-- Reference: feature_sets_studies (table: feature_sets)
ALTER TABLE feature_sets ADD CONSTRAINT feature_sets_studies FOREIGN KEY feature_sets_studies (studyID)
    REFERENCES studies (studyID);

-- Reference: linear_model_session_lm_analysis (table: linear_model)
ALTER TABLE linear_model ADD CONSTRAINT linear_model_session_lm_analysis FOREIGN KEY linear_model_session_lm_analysis (sessionID)
    REFERENCES session_lm_analysis (session_analysis_ID);

-- Reference: lm_cohend_results_lm_results_keys (table: lm_cohend_results)
ALTER TABLE lm_cohend_results ADD CONSTRAINT lm_cohend_results_lm_results_keys FOREIGN KEY lm_cohend_results_lm_results_keys (res_keyID)
    REFERENCES lm_results_keys (res_keyID);

-- Reference: lm_corr_results_lm_results_keys (table: lm_corr_results)
ALTER TABLE lm_corr_results ADD CONSTRAINT lm_corr_results_lm_results_keys FOREIGN KEY lm_corr_results_lm_results_keys (res_keyID)
    REFERENCES lm_results_keys (res_keyID);

-- Reference: lm_demog_results_lm_results_keys (table: lm_demog_results)
ALTER TABLE lm_demog_results ADD CONSTRAINT lm_demog_results_lm_results_keys FOREIGN KEY lm_demog_results_lm_results_keys (res_keyID)
    REFERENCES lm_results_keys (res_keyID);

-- Reference: lm_filter_linear_model (table: lm_filter)
ALTER TABLE lm_filter ADD CONSTRAINT lm_filter_linear_model FOREIGN KEY lm_filter_linear_model (lmID,sessionID)
    REFERENCES linear_model (lmID,sessionID);

-- Reference: lm_interactions_linear_model (table: lm_interactions)
ALTER TABLE lm_interactions ADD CONSTRAINT lm_interactions_linear_model FOREIGN KEY lm_interactions_linear_model (lmID,sessionID)
    REFERENCES linear_model (lmID,sessionID);

-- Reference: lm_mutate_linear_model (table: lm_mutate)
ALTER TABLE lm_mutate ADD CONSTRAINT lm_mutate_linear_model FOREIGN KEY lm_mutate_linear_model (lmID,sessionID)
    REFERENCES linear_model (lmID,sessionID);

-- Reference: lm_results_keys_linear_model (table: lm_results_keys)
ALTER TABLE lm_results_keys ADD CONSTRAINT lm_results_keys_linear_model FOREIGN KEY lm_results_keys_linear_model (lmID,sessionID)
    REFERENCES linear_model (lmID,sessionID);

-- Reference: lm_results_lm_results_keys (table: lm_results)
ALTER TABLE lm_results ADD CONSTRAINT lm_results_lm_results_keys FOREIGN KEY lm_results_lm_results_keys (res_keyID)
    REFERENCES lm_results_keys (res_keyID);

-- Reference: lm_summary_fac_lm_results_keys (table: lm_summary_fact)
ALTER TABLE lm_summary_fact ADD CONSTRAINT lm_summary_fac_lm_results_keys FOREIGN KEY lm_summary_fac_lm_results_keys (res_keyID)
    REFERENCES lm_results_keys (res_keyID);

-- Reference: lm_summary_lm_results_keys (table: lm_summary_cont)
ALTER TABLE lm_summary_cont ADD CONSTRAINT lm_summary_lm_results_keys FOREIGN KEY lm_summary_lm_results_keys (res_keyID)
    REFERENCES lm_results_keys (res_keyID);

-- Reference: lm_variables_linear_model (table: lm_variables)
ALTER TABLE lm_variables ADD CONSTRAINT lm_variables_linear_model FOREIGN KEY lm_variables_linear_model (lmID,sessionID)
    REFERENCES linear_model (lmID,sessionID);

-- Reference: meta_beta_results_session_meta (table: meta_beta_results)
ALTER TABLE meta_beta_results ADD CONSTRAINT meta_beta_results_session_meta FOREIGN KEY meta_beta_results_session_meta (meta_key_ID)
    REFERENCES session_meta (meta_key_ID);

-- Reference: meta_cohd_results_session_meta (table: meta_cohd_results)
ALTER TABLE meta_cohd_results ADD CONSTRAINT meta_cohd_results_session_meta FOREIGN KEY meta_cohd_results_session_meta (meta_key_ID)
    REFERENCES session_meta (meta_key_ID);

-- Reference: metrics_data_session_loadMetrics (table: metrics_data)
ALTER TABLE metrics_data ADD CONSTRAINT metrics_data_session_loadMetrics FOREIGN KEY metrics_data_session_loadMetrics (session_metr_ID)
    REFERENCES session_loadMetrics (session_metr_ID);

-- Reference: metrics_timestamp_session_loadMetrics (table: metrics_file)
ALTER TABLE metrics_file ADD CONSTRAINT metrics_timestamp_session_loadMetrics FOREIGN KEY metrics_timestamp_session_loadMetrics (session_metr_ID)
    REFERENCES session_loadMetrics (session_metr_ID);

-- Reference: metrics_timestamp_sites_in_study (table: metrics_file)
ALTER TABLE metrics_file ADD CONSTRAINT metrics_timestamp_sites_in_study FOREIGN KEY metrics_timestamp_sites_in_study (siteID,studyID)
    REFERENCES sites_in_study (siteID,studyID);

-- Reference: metrics_timestamp_studies (table: metrics_file)
ALTER TABLE metrics_file ADD CONSTRAINT metrics_timestamp_studies FOREIGN KEY metrics_timestamp_studies (studyID)
    REFERENCES studies (studyID);

-- Reference: session_covariates_sites_in_study (table: session_covariates)
ALTER TABLE session_covariates ADD CONSTRAINT session_covariates_sites_in_study FOREIGN KEY session_covariates_sites_in_study (study_site_ID)
    REFERENCES sites_in_study (study_site_ID);

-- Reference: session_lm_analysis_studies (table: session_lm_analysis)
ALTER TABLE session_lm_analysis ADD CONSTRAINT session_lm_analysis_studies FOREIGN KEY session_lm_analysis_studies (studyID)
    REFERENCES studies (studyID);

-- Reference: session_loadMetrics_sites_in_study (table: session_loadMetrics)
ALTER TABLE session_loadMetrics ADD CONSTRAINT session_loadMetrics_sites_in_study FOREIGN KEY session_loadMetrics_sites_in_study (study_site_ID)
    REFERENCES sites_in_study (study_site_ID);

-- Reference: session_meta_studies (table: session_meta)
ALTER TABLE session_meta ADD CONSTRAINT session_meta_studies FOREIGN KEY session_meta_studies (studyID)
    REFERENCES studies (studyID);

-- Reference: sites_in_meta_session_meta (table: sites_in_meta)
ALTER TABLE sites_in_meta ADD CONSTRAINT sites_in_meta_session_meta FOREIGN KEY sites_in_meta_session_meta (session_meta_ID,lmID,ROI)
    REFERENCES session_meta (session_meta_ID,lmID,ROI);

-- Reference: sites_studies (table: sites_in_study)
ALTER TABLE sites_in_study ADD CONSTRAINT sites_studies FOREIGN KEY sites_studies (studyID)
    REFERENCES studies (studyID);

-- Reference: subj_excluded_session_loadMetrics (table: subj_excluded)
ALTER TABLE subj_excluded ADD CONSTRAINT subj_excluded_session_loadMetrics FOREIGN KEY subj_excluded_session_loadMetrics (session_metr_ID)
    REFERENCES session_loadMetrics (session_metr_ID);

-- Reference: subj_excluded_sites_in_study (table: subj_excluded)
ALTER TABLE subj_excluded ADD CONSTRAINT subj_excluded_sites_in_study FOREIGN KEY subj_excluded_sites_in_study (siteID,studyID)
    REFERENCES sites_in_study (siteID,studyID);

-- End of file.

