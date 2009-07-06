/*
 * parseInv.cpp
 *
 *  Created on: 18 juin 2009
 *      Author: matthieu
 */

#include "parseInv.h"

int parse_saneInv_ini_file(char * ini_name, string &fname,string &boloname, string &noiseSp_dir_output, string &extentnoiseSp)
{
	dictionary	*	ini ;

	/* Some temporary variables to hold query results */
	//int				b ;
	//int				i ;
	//double			d ;
	char		*	s ;
	//long l;
	//string str;

	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	s = iniparser_getstring(ini, "commons:channel",NULL);
	if(s!=NULL){
		printf("channel file : [%s]\n",s);
		boloname=s;
		//read_bolofile((string)s, bolonames);
	}else{
		printf("You must specify a bolometer file : commons:channel\n");
		return(-1);
	}//	channel =./RCW_120_M/bolos_commons.txt ;

	s = iniparser_getstring(ini, "sanepic_inv_matrix:cov_matrix_file",NULL);
	if(s!=NULL){
		printf("cov_matrix_file: [%s]\n",s);
		fname=s;
		//read_bolofile((string)s, bolonames);
	}else{
		printf("You must specify a noise covariance matrix file to invert : sanepic_inv_matrix:cov_matrix_file\n");
		return(-1);
	}//	fname = ./RCW_120_M/BoloPS0sanepic_binary.psd

	s = iniparser_getstring(ini, "sanepic_inv_matrix:noise_out_dir",NULL);
	if(s!=NULL){
		printf("noise_out_dir: [%s]\n",s);
		noiseSp_dir_output=s;
		//read_bolofile((string)s, bolonames);
	}else{
		printf("You must specify a output directory : sanepic_inv_matrix:noise_out_dir\n");
		return(-1);
	}//	noiseSp_dir_output = ./RCW_120_M/

	s = iniparser_getstring(ini, "sanepic_inv_matrix:noise_prefixe",NULL);
	if(s!=NULL){
		printf("noise_prefixe: [%s]\n",s);
		extentnoiseSp=s;
		//read_bolofile((string)s, bolonames);
	}else{
		printf("You must specify a output directory : sanepic_inv_matrix:noise_prefixe\n");
		return(-1);
	}//	extentnoiseSp = NoisePS


	iniparser_freedict(ini);

	return 0;
}