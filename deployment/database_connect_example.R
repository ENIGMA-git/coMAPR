database_connect <- function () {
	con <- eregr_connect("db_name_obsolete",RMySQL::MySQL(),user='username',password='password',dbname='enigma_par_db_name',host='localhost or host name',port=3306)
	con
}

