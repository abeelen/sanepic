/*
 * parseInv.cpp
 *
 *  Created on: 18 juin 2009
 *      Author: matthieu
 */

#include "parseInv.h"


//#include <tchar.h>
//#include <windows.h>
//#include <stdlib.h>
//#include <stdio.h>
#include <iostream>
#include <string>

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}


using namespace std;

std::string		Basename(std::string path);

int parse_saneInv_ini_file(char * ini_name, string &fname,string &boloname, string &noiseSp_dir_output, string &extentnoiseSp)
{
	dictionary	*	ini ;

	/* Some temporary variables to hold query results */
	char		*	s ;


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

	s = iniparser_getstring(ini, "sanepic_inv_matrix:cov_matrix_file",NULL);
	if(s!=NULL){
		printf("cov_matrix_file: [%s]\n",s);
		fname=s;
//		Basename(fname);
		//read_strings((string)s, bolonames);
	}else{
		printf("You must specify a noise covariance matrix file to invert : sanepic_inv_matrix:cov_matrix_file\n");
		return(-1);
	}//	fname = ./RCW_120_M/BoloPS0sanepic_binary.psd

	char * pPath;
	pPath = getenv ("TMPBATCH");
	if (pPath!=NULL){
		noiseSp_dir_output=pPath;
		printf ("The current path is: %s\n",pPath);
	}else{
		s = iniparser_getstring(ini, "commons:temp_dir",NULL);
		if(s!=NULL){
			printf("noise_out_dir: [%s]\n",s);
			noiseSp_dir_output=s;
			//read_strings((string)s, bolonames);
		}else{
			printf("You must specify a output directory : commons:temp_dir\n");
			return(-1);
		}//	noiseSp_dir_output = ./RCW_120_M/
	}


	s = iniparser_getstring(ini, "sanepic_inv_matrix:noise_prefixe",NULL); // TODO : change to noise_prefixe_file
	if(s!=NULL){
		printf("noise_prefixe: [%s]\n",s);
		extentnoiseSp=s;
		//read_strings((string)s, bolonames);
	}else{
		printf("You must specify a output directory : sanepic_inv_matrix:noise_prefixe\n");
		return(-1);
	}//	extentnoiseSp = NoisePS

	// cleaning up
	iniparser_freedict(ini);

	return 0;
}


std::string		Basename(std::string path)
{
	std::string		Result;



	char * pch;
	char * temp=NULL;
	printf ("Splitting string \"%s\" into tokens:\n",path.c_str());
	pch = strtok ((char*)path.c_str(),".");
	while (pch != NULL)
	{
		if (strcmp (pch,(char*)"fits") == 0)
			break;
		temp=pch;
		printf ("%s\n",temp);
		pch = strtok (NULL, ".");

	}
	cout << "result : " << temp << endl;
	getchar();
	//	Result = fname;
	//	Result += ext;
	return Result;

}

