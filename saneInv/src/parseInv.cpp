/*
 * parseInv.cpp
 *
 *  Created on: 18 juin 2009
 *      Author: matthieu
 */

#include "parseInv.h"
#include "inputFileIO.h"

//#include <tchar.h>
//#include <windows.h>
//#include <stdlib.h>
//#include <stdio.h>
#include "parser_functions.h"
#include <iostream>
#include <string>

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

using namespace std;


int parse_saneInv_ini_file(char * ini_name, string &fname,struct samples &samples_struct,struct directories &dir,string &boloname, string &base)
{
	dictionary	*	ini ;

	/* Some temporary variables to hold query results */
	char		*	s ;
	//	string base;


	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) { // if dictionnary was not found, exit + error msg
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	s = iniparser_getstring(ini, "commons:channel",NULL);
	if(s!=NULL){
		//		printf("channel file : [%s]\n",s);
		boloname=s;
		//read_strings((string)s, bolonames);
	}else{
		printf("You must specify a bolometer file : commons:channel\n");
		return(-1);
	}//	channel =./RCW_120_M/bolos_commons.txt ;

	s = iniparser_getstring(ini, "commons:noise_prefixe",NULL);
	if(s!=NULL){
		fname=s;
		//		base=Basename(fname);
	}else{
		printf("Warning ! You must specify a noise covariance matrix file to invert : commons:noise_prefixe or put the information in the fits filelist\n");
		//		return(-1);
	}//	fname = ./RCW_120_M/BoloPS0sanepic_binary.psd


	if(read_directories(ini, dir, 0)==-1)
		return -1;


	if(read_fits_file_list(ini, dir,samples_struct, 0)==-1)
		return -1;

//	printf("\nsaneInv parser operations completed :\n");
	cout << "You have specified the following options : \n";

	print_directories(dir);
	//printf("cov_matrix_file: [%s]\n",s);

	// cleaning up
	iniparser_freedict(ini);

	return 0;
}


