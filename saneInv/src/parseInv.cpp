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
		printf("channel file : [%s]\n",s);
		boloname=s;
		//read_strings((string)s, bolonames);
	}else{
		printf("You must specify a bolometer file : commons:channel\n");
		return(-1);
	}//	channel =./RCW_120_M/bolos_commons.txt ;

//	s = iniparser_getstring(ini, "sanepic_inv_matrix:cov_matrix_file",NULL);
	s = iniparser_getstring(ini, "commons:noise_prefixe",NULL);
	if(s!=NULL){
		printf("cov_matrix_file: [%s]\n",s);
		fname=s;
//		base=Basename(fname);
	}else{
		printf("Warning ! You must specify a noise covariance matrix file to invert : commons:noise_prefixe or put the information in the fits filelist\n");
//		return(-1);
	}//	fname = ./RCW_120_M/BoloPS0sanepic_binary.psd

//	char * pPath;
//	pPath = getenv ("TMPBATCH");
//	if (pPath!=NULL){
//		noiseSp_dir_output=pPath;
//		printf ("The current path is: %s\n",pPath);
//	}else{
//		s = iniparser_getstring(ini, "commons:temp_dir",NULL);
//		if(s!=NULL){
//			printf("noise_out_dir: [%s]\n",s);
//			noiseSp_dir_output=s;
//			//read_strings((string)s, bolonames);
//		}else{
//			printf("You must specify a output directory : commons:temp_dir\n");
//			return(-1);
//		}//	noiseSp_dir_output = ./RCW_120_M/
//	}


	if(read_directories(ini, dir, 0)==-1)
		return -1;


	if(read_fits_file_list(ini, dir,samples_struct, 0)==-1)
		return -1;


	//	s = iniparser_getstring(ini, "sanepic_inv_matrix:noise_prefixe",NULL);
	//	if(s!=NULL){
	//		printf("noise_prefixe: [%s]\n",s);
	//		extentnoiseSp=s;
	//		//read_strings((string)s, bolonames);
	//	}else{
	//		printf("You must specify a output directory : sanepic_inv_matrix:noise_prefixe\n");
	//		return(-1);
	//	}//	extentnoiseSp = NoisePS

	// cleaning up
	iniparser_freedict(ini);

	return 0;
}


