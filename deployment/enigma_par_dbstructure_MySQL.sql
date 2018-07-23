-- MySQL dump 10.13  Distrib 5.7.17, for macos10.12 (x86_64)
--
-- Host: igc-mysql.ini.usc.edu    Database: enigma_mdddti
-- ------------------------------------------------------
-- Server version	5.7.18

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `covariates_general`
--

DROP TABLE IF EXISTS `covariates_general`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `covariates_general` (
  `subjID` char(100) NOT NULL,
  `cov_name` char(100) NOT NULL,
  `cov_value` char(255) DEFAULT NULL,
  `session_covar_ID` char(200) NOT NULL,
  PRIMARY KEY (`subjID`,`cov_name`,`session_covar_ID`),
  KEY `covariates_general_idx_1` (`subjID`,`cov_name`),
  KEY `covariates_general_session_covariates` (`session_covar_ID`),
  CONSTRAINT `covariates_general_session_covariates` FOREIGN KEY (`session_covar_ID`) REFERENCES `session_covariates` (`session_covar_ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `feature_sets`
--

DROP TABLE IF EXISTS `feature_sets`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `feature_sets` (
  `session_fs_ID` char(200) NOT NULL,
  `fsID` char(100) NOT NULL,
  `metric` char(100) NOT NULL,
  `ROI` char(100) NOT NULL,
  PRIMARY KEY (`session_fs_ID`,`fsID`,`metric`,`ROI`),
  CONSTRAINT `feature_sets_session_fs` FOREIGN KEY (`session_fs_ID`, `fsID`) REFERENCES `session_fs` (`session_fs_ID`, `fsID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `feature_sets_covariates`
--

DROP TABLE IF EXISTS `feature_sets_covariates`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `feature_sets_covariates` (
  `session_fs_ID` char(200) NOT NULL,
  `fsID` char(100) NOT NULL,
  `subjID` char(100) NOT NULL,
  `cov_name` char(100) NOT NULL,
  `cov_value` double NOT NULL,
  `metric` char(100) NOT NULL,
  `vertex` int(11) NOT NULL,
  `row_ID` int(11) NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (`session_fs_ID`,`fsID`,`row_ID`),
  UNIQUE KEY `row_ID_KEY` (`row_ID`),
  CONSTRAINT `feature_sets_covariates_session_fs` FOREIGN KEY (`session_fs_ID`, `fsID`) REFERENCES `session_fs` (`session_fs_ID`, `fsID`)
) ENGINE=InnoDB AUTO_INCREMENT=316674 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `feature_sets_newregr`
--

DROP TABLE IF EXISTS `feature_sets_newregr`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `feature_sets_newregr` (
  `session_fs_ID` char(200) NOT NULL,
  `fsID` char(100) NOT NULL,
  `var` char(100) NOT NULL,
  `formula` varchar(1000) NOT NULL,
  PRIMARY KEY (`session_fs_ID`,`fsID`,`var`),
  CONSTRAINT `Table_30_session_fs` FOREIGN KEY (`session_fs_ID`, `fsID`) REFERENCES `session_fs` (`session_fs_ID`, `fsID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `linear_model`
--

DROP TABLE IF EXISTS `linear_model`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `linear_model` (
  `lmID` char(100) NOT NULL,
  `sessionID` char(200) NOT NULL,
  `fsID` char(100) NOT NULL,
  `session_fs_ID` char(200) NOT NULL,
  `lm_gDoc_path` varchar(1000) NOT NULL,
  `lm_name` varchar(1000) NOT NULL,
  `lm_text` varchar(1000) NOT NULL,
  `main_factor` char(100) DEFAULT NULL,
  `new_regressors` varchar(1000) NOT NULL,
  `cont_value` int(11) DEFAULT '0',
  `pat_value` int(11) DEFAULT '1',
  `cont_min` int(11) DEFAULT '0',
  `pat_min` int(11) DEFAULT '0',
  `comments` varchar(1000) DEFAULT NULL,
  PRIMARY KEY (`lmID`,`sessionID`),
  KEY `linear_model_idx_1` (`sessionID`),
  KEY `linear_model_idx_2` (`lmID`,`sessionID`),
  KEY `linear_model_session_fs` (`session_fs_ID`,`fsID`),
  CONSTRAINT `linear_model_session_fs` FOREIGN KEY (`session_fs_ID`, `fsID`) REFERENCES `session_fs` (`session_fs_ID`, `fsID`),
  CONSTRAINT `linear_model_session_lm_analysis` FOREIGN KEY (`sessionID`) REFERENCES `session_lm_analysis` (`session_analysis_ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lm_cohend_results`
--

DROP TABLE IF EXISTS `lm_cohend_results`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `lm_cohend_results` (
  `res_keyID` bigint(20) NOT NULL,
  `vertex` int(11) NOT NULL,
  `var` char(100) NOT NULL,
  `cohens_d` double NOT NULL,
  `cohens_se` double NOT NULL,
  `cohens_low_ci` double NOT NULL,
  `cohens_high_ci` double NOT NULL,
  `cohens_pval` double NOT NULL,
  PRIMARY KEY (`res_keyID`,`vertex`,`var`),
  CONSTRAINT `lm_cohend_results_lm_results_keys` FOREIGN KEY (`res_keyID`) REFERENCES `lm_results_keys` (`res_keyID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lm_corr_results`
--

DROP TABLE IF EXISTS `lm_corr_results`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `lm_corr_results` (
  `res_keyID` bigint(20) NOT NULL,
  `vertex` char(100) NOT NULL,
  `var` char(100) NOT NULL,
  `corr` double NOT NULL,
  `corr_pval` double NOT NULL,
  `corr_se` double NOT NULL,
  PRIMARY KEY (`res_keyID`,`vertex`,`var`),
  CONSTRAINT `lm_corr_results_lm_results_keys` FOREIGN KEY (`res_keyID`) REFERENCES `lm_results_keys` (`res_keyID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lm_demog_results`
--

DROP TABLE IF EXISTS `lm_demog_results`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `lm_demog_results` (
  `res_keyID` bigint(20) NOT NULL,
  `n_overall` int(11) NOT NULL,
  `n_cont` int(11) NOT NULL,
  `n_pat` int(11) NOT NULL,
  PRIMARY KEY (`res_keyID`),
  CONSTRAINT `lm_demog_results_lm_results_keys` FOREIGN KEY (`res_keyID`) REFERENCES `lm_results_keys` (`res_keyID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lm_filter`
--

DROP TABLE IF EXISTS `lm_filter`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `lm_filter` (
  `lmID` char(100) NOT NULL,
  `sessionID` char(200) NOT NULL,
  `filter_1` varchar(1000) DEFAULT NULL,
  `filter_2` varchar(1000) DEFAULT NULL,
  `filter_full` varchar(1000) DEFAULT NULL,
  PRIMARY KEY (`lmID`,`sessionID`),
  KEY `lm_filter_idx_1` (`lmID`,`sessionID`),
  CONSTRAINT `lm_filter_linear_model` FOREIGN KEY (`lmID`, `sessionID`) REFERENCES `linear_model` (`lmID`, `sessionID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lm_interactions`
--

DROP TABLE IF EXISTS `lm_interactions`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `lm_interactions` (
  `var1` char(100) NOT NULL,
  `var2` char(100) NOT NULL,
  `lmID` char(100) NOT NULL,
  `sessionID` char(100) NOT NULL,
  PRIMARY KEY (`var1`,`var2`,`lmID`,`sessionID`),
  KEY `lm_interactions_idx_1` (`lmID`,`sessionID`),
  CONSTRAINT `lm_interactions_linear_model` FOREIGN KEY (`lmID`, `sessionID`) REFERENCES `linear_model` (`lmID`, `sessionID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lm_mutate`
--

DROP TABLE IF EXISTS `lm_mutate`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `lm_mutate` (
  `formula` varchar(1000) NOT NULL,
  `order` int(11) NOT NULL,
  `var` char(100) NOT NULL,
  `lmID` char(100) NOT NULL,
  `sessionID` char(200) NOT NULL,
  PRIMARY KEY (`var`,`lmID`,`sessionID`),
  KEY `lm_mutate_idx_1` (`lmID`,`sessionID`),
  CONSTRAINT `lm_mutate_linear_model` FOREIGN KEY (`lmID`, `sessionID`) REFERENCES `linear_model` (`lmID`, `sessionID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lm_results`
--

DROP TABLE IF EXISTS `lm_results`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `lm_results` (
  `res_keyID` bigint(20) NOT NULL,
  `vertex` int(11) NOT NULL,
  `var` char(100) NOT NULL,
  `beta` double NOT NULL,
  `sterr` double NOT NULL,
  `p_beta` double NOT NULL,
  PRIMARY KEY (`res_keyID`,`vertex`,`var`),
  CONSTRAINT `lm_results_lm_results_keys` FOREIGN KEY (`res_keyID`) REFERENCES `lm_results_keys` (`res_keyID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lm_results_keys`
--

DROP TABLE IF EXISTS `lm_results_keys`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `lm_results_keys` (
  `res_keyID` bigint(20) NOT NULL AUTO_INCREMENT,
  `lmID` char(100) NOT NULL,
  `sessionID` char(200) NOT NULL,
  `metric` char(100) NOT NULL,
  `ROI` char(100) NOT NULL,
  `result_sessionID` char(200) NOT NULL,
  PRIMARY KEY (`res_keyID`),
  KEY `lm_results_keys_linear_model` (`lmID`,`sessionID`),
  CONSTRAINT `lm_results_keys_linear_model` FOREIGN KEY (`lmID`, `sessionID`) REFERENCES `linear_model` (`lmID`, `sessionID`)
) ENGINE=InnoDB AUTO_INCREMENT=1333108 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lm_summary_cont`
--

DROP TABLE IF EXISTS `lm_summary_cont`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `lm_summary_cont` (
  `res_keyID` bigint(20) NOT NULL,
  `var` char(100) NOT NULL,
  `mean` double NOT NULL,
  `median` double NOT NULL,
  `Q1` double NOT NULL,
  `Q3` double NOT NULL,
  `min` double NOT NULL,
  `max` double NOT NULL,
  `std` double NOT NULL,
  PRIMARY KEY (`res_keyID`,`var`),
  CONSTRAINT `lm_summary_lm_results_keys` FOREIGN KEY (`res_keyID`) REFERENCES `lm_results_keys` (`res_keyID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lm_summary_fact`
--

DROP TABLE IF EXISTS `lm_summary_fact`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `lm_summary_fact` (
  `res_keyID` bigint(20) NOT NULL,
  `var` char(100) NOT NULL,
  `value` int(11) NOT NULL,
  `amount` int(11) NOT NULL,
  PRIMARY KEY (`res_keyID`,`var`,`value`),
  CONSTRAINT `lm_summary_fac_lm_results_keys` FOREIGN KEY (`res_keyID`) REFERENCES `lm_results_keys` (`res_keyID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lm_variables`
--

DROP TABLE IF EXISTS `lm_variables`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `lm_variables` (
  `var` char(100) NOT NULL,
  `lmID` char(100) NOT NULL,
  `sessionID` char(200) NOT NULL,
  `modifier` char(100) DEFAULT NULL,
  `is_global` tinyint(1) NOT NULL,
  PRIMARY KEY (`var`,`lmID`,`sessionID`),
  KEY `lm_variables_idx_1` (`lmID`,`sessionID`),
  CONSTRAINT `lm_variables_linear_model` FOREIGN KEY (`lmID`, `sessionID`) REFERENCES `linear_model` (`lmID`, `sessionID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `meta_beta_results`
--

DROP TABLE IF EXISTS `meta_beta_results`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `meta_beta_results` (
  `meta_key_ID` int(11) NOT NULL,
  `metric` char(100) NOT NULL,
  `vertex` int(11) NOT NULL,
  `var` char(100) NOT NULL,
  `b` double NOT NULL,
  `se` double NOT NULL,
  `zval` double NOT NULL,
  `pval` double NOT NULL,
  `ci_lb` double NOT NULL,
  `ci_ub` double NOT NULL,
  `tau2` double NOT NULL,
  `se_tau2` double NOT NULL,
  `I2` double NOT NULL,
  `H2` double NOT NULL,
  PRIMARY KEY (`meta_key_ID`,`metric`,`vertex`,`var`),
  CONSTRAINT `meta_beta_results_session_meta` FOREIGN KEY (`meta_key_ID`) REFERENCES `session_meta` (`meta_key_ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `meta_cohd_results`
--

DROP TABLE IF EXISTS `meta_cohd_results`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `meta_cohd_results` (
  `meta_key_ID` int(11) NOT NULL,
  `metric` char(100) NOT NULL,
  `vertex` int(11) NOT NULL,
  `var` char(100) NOT NULL,
  `cohd` double NOT NULL,
  `se` double NOT NULL,
  `zval` double NOT NULL,
  `pval` double NOT NULL,
  `ci_lb` double NOT NULL,
  `ci_ub` double NOT NULL,
  `tau2` double NOT NULL,
  `se_tau2` double NOT NULL,
  `I2` double NOT NULL,
  `H2` double NOT NULL,
  PRIMARY KEY (`meta_key_ID`,`metric`,`vertex`,`var`),
  CONSTRAINT `meta_cohd_results_session_meta` FOREIGN KEY (`meta_key_ID`) REFERENCES `session_meta` (`meta_key_ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `metrics_data`
--

DROP TABLE IF EXISTS `metrics_data`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `metrics_data` (
  `session_metr_ID` char(200) NOT NULL,
  `subjID` char(100) NOT NULL,
  `metric` char(100) NOT NULL,
  `ROI` char(100) NOT NULL,
  `vertex` int(11) NOT NULL,
  `value` double DEFAULT NULL,
  `row_ID` int(11) NOT NULL AUTO_INCREMENT,
  `fsID` char(100) DEFAULT NULL,
  PRIMARY KEY (`session_metr_ID`,`row_ID`),
  UNIQUE KEY `metrics_data_row_ID_KEY` (`row_ID`),
  KEY `metrics_data_idx_1` (`subjID`,`metric`,`ROI`),
  KEY `metrics_data_idx_2` (`row_ID`),
  CONSTRAINT `metrics_data_session_loadMetrics` FOREIGN KEY (`session_metr_ID`) REFERENCES `session_loadMetrics` (`session_metr_ID`)
) ENGINE=InnoDB AUTO_INCREMENT=35052186 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `metrics_file`
--

DROP TABLE IF EXISTS `metrics_file`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `metrics_file` (
  `studyID` char(100) NOT NULL,
  `siteID` char(100) NOT NULL,
  `subjID` char(100) NOT NULL,
  `metric` char(100) NOT NULL,
  `file_path` varchar(1000) NOT NULL,
  `file_lastedit` datetime NOT NULL,
  `ROI` char(100) NOT NULL,
  `session_metr_ID` char(200) NOT NULL,
  PRIMARY KEY (`studyID`,`siteID`,`subjID`,`metric`,`ROI`),
  KEY `metrics_file_idx_1` (`studyID`,`siteID`,`subjID`,`ROI`),
  KEY `metrics_timestamp_session_loadMetrics` (`session_metr_ID`),
  KEY `metrics_timestamp_sites_in_study` (`siteID`,`studyID`),
  CONSTRAINT `metrics_timestamp_session_loadMetrics` FOREIGN KEY (`session_metr_ID`) REFERENCES `session_loadMetrics` (`session_metr_ID`),
  CONSTRAINT `metrics_timestamp_sites_in_study` FOREIGN KEY (`siteID`, `studyID`) REFERENCES `sites_in_study` (`siteID`, `studyID`),
  CONSTRAINT `metrics_timestamp_studies` FOREIGN KEY (`studyID`) REFERENCES `studies` (`studyID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `session_covariates`
--

DROP TABLE IF EXISTS `session_covariates`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `session_covariates` (
  `session_covar_ID` char(200) NOT NULL,
  `timestamp` datetime NOT NULL,
  `file_cov_path` varchar(1000) NOT NULL,
  `file_lastedit` datetime NOT NULL,
  `study_site_ID` int(11) NOT NULL,
  PRIMARY KEY (`session_covar_ID`),
  KEY `session_covariates_idx_1` (`session_covar_ID`),
  KEY `session_covariates_sites_in_study` (`study_site_ID`),
  CONSTRAINT `session_covariates_sites_in_study` FOREIGN KEY (`study_site_ID`) REFERENCES `sites_in_study` (`study_site_ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `session_fs`
--

DROP TABLE IF EXISTS `session_fs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `session_fs` (
  `session_fs_ID` char(200) NOT NULL,
  `fsID` char(100) NOT NULL,
  `timestamp` datetime NOT NULL,
  `study_site_ID` int(11) NOT NULL,
  PRIMARY KEY (`session_fs_ID`,`fsID`),
  KEY `session_fs_sites_in_study` (`study_site_ID`),
  CONSTRAINT `session_fs_sites_in_study` FOREIGN KEY (`study_site_ID`) REFERENCES `sites_in_study` (`study_site_ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `session_lm_analysis`
--

DROP TABLE IF EXISTS `session_lm_analysis`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `session_lm_analysis` (
  `session_analysis_ID` char(200) NOT NULL,
  `studyID` char(100) NOT NULL,
  `siteID` char(100) NOT NULL,
  `timestamp` datetime NOT NULL,
  PRIMARY KEY (`session_analysis_ID`),
  KEY `session_lm_analysis_idx_1` (`session_analysis_ID`),
  KEY `session_lm_analysis_studies` (`studyID`),
  CONSTRAINT `session_lm_analysis_studies` FOREIGN KEY (`studyID`) REFERENCES `studies` (`studyID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `session_loadMetrics`
--

DROP TABLE IF EXISTS `session_loadMetrics`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `session_loadMetrics` (
  `session_metr_ID` char(200) NOT NULL,
  `timestamp` datetime NOT NULL,
  `study_site_ID` int(11) NOT NULL,
  PRIMARY KEY (`session_metr_ID`),
  KEY `session_loadMetrics_sites_in_study` (`study_site_ID`),
  CONSTRAINT `session_loadMetrics_sites_in_study` FOREIGN KEY (`study_site_ID`) REFERENCES `sites_in_study` (`study_site_ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `session_meta`
--

DROP TABLE IF EXISTS `session_meta`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `session_meta` (
  `meta_key_ID` int(11) NOT NULL AUTO_INCREMENT,
  `session_meta_ID` char(200) NOT NULL,
  `studyID` char(100) NOT NULL,
  `lmID` char(100) NOT NULL,
  `ROI` char(100) NOT NULL,
  PRIMARY KEY (`meta_key_ID`),
  UNIQUE KEY `session_meta_ak_1` (`session_meta_ID`,`lmID`,`ROI`),
  KEY `session_meta_studies` (`studyID`),
  CONSTRAINT `session_meta_studies` FOREIGN KEY (`studyID`) REFERENCES `studies` (`studyID`)
) ENGINE=InnoDB AUTO_INCREMENT=7034 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `sites_in_meta`
--

DROP TABLE IF EXISTS `sites_in_meta`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sites_in_meta` (
  `session_meta_ID` char(200) NOT NULL,
  `lmID` char(100) NOT NULL,
  `ROI` char(100) NOT NULL,
  `siteID` char(100) NOT NULL,
  `result_sessionID` char(200) NOT NULL,
  PRIMARY KEY (`session_meta_ID`,`lmID`,`ROI`,`siteID`),
  CONSTRAINT `sites_in_meta_session_meta` FOREIGN KEY (`session_meta_ID`, `lmID`, `ROI`) REFERENCES `session_meta` (`session_meta_ID`, `lmID`, `ROI`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `sites_in_study`
--

DROP TABLE IF EXISTS `sites_in_study`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sites_in_study` (
  `siteID` char(100) NOT NULL,
  `studyID` char(100) NOT NULL,
  `study_site_ID` int(11) NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (`siteID`,`studyID`),
  UNIQUE KEY `study_siteID` (`study_site_ID`),
  KEY `sites_studies` (`studyID`),
  CONSTRAINT `sites_studies` FOREIGN KEY (`studyID`) REFERENCES `studies` (`studyID`)
) ENGINE=InnoDB AUTO_INCREMENT=626 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `studies`
--

DROP TABLE IF EXISTS `studies`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `studies` (
  `studyID` char(100) NOT NULL,
  `study_name` varchar(100) DEFAULT NULL,
  `study_gDoc_path` varchar(1000) NOT NULL,
  `study_analysis_path` varchar(1000) NOT NULL,
  `study_fs_path` varchar(1000) NOT NULL,
  `study_demographics_path` varchar(1000) NOT NULL,
  `study_data_format` char(100) NOT NULL,
  PRIMARY KEY (`studyID`),
  KEY `studies_idx_1` (`studyID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `study_metrics`
--

DROP TABLE IF EXISTS `study_metrics`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `study_metrics` (
  `studyID` char(100) NOT NULL,
  `metr_name` char(100) NOT NULL,
  PRIMARY KEY (`studyID`,`metr_name`),
  KEY `study_metrics_idx_1` (`studyID`,`metr_name`),
  CONSTRAINT `Table_14_studies` FOREIGN KEY (`studyID`) REFERENCES `studies` (`studyID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `subj_excluded`
--

DROP TABLE IF EXISTS `subj_excluded`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `subj_excluded` (
  `siteID` char(100) NOT NULL,
  `studyID` char(100) NOT NULL,
  `subjID` char(100) NOT NULL,
  `ROI` char(100) NOT NULL,
  `QC` int(11) DEFAULT NULL,
  `session_metr_ID` char(200) NOT NULL,
  PRIMARY KEY (`siteID`,`studyID`,`subjID`,`ROI`),
  KEY `subj_excluded_session_loadMetrics` (`session_metr_ID`),
  CONSTRAINT `subj_excluded_session_loadMetrics` FOREIGN KEY (`session_metr_ID`) REFERENCES `session_loadMetrics` (`session_metr_ID`),
  CONSTRAINT `subj_excluded_sites_in_study` FOREIGN KEY (`siteID`, `studyID`) REFERENCES `sites_in_study` (`siteID`, `studyID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2018-07-23 11:45:40
