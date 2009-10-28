/*
 * parsePS.cpp
 *
 *  Created on: 18 juin 2009
 *      Author: matthieu
 */


#include <iostream>
#include <string>
#include <vector>

#include "parsePS.h"
#include "mpi_architecture_builder.h"
#include "parser_functions.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

using namespace std;

int parse_sanePS_ini_file(char * ini_name, struct user_options &u_opt, struct directories &dir, struct samples &samples_struct,struct input_commons &com,
		struct detectors &det,string &MixMatfile, string &ellFile, string &signame)

{
	dictionary	*	ini ;

	//int nnf=0; // extentnoiseSp_list number of elements


	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	// printf dictionnary to stderr for debugging
	//iniparser_dump(ini, stderr);


	/* Get sanepic_preprocess attributes */
	printf("\nSanepic Noise Estimation Procedure:\n");


	if(read_directories(ini, dir)==-1)
		return -1;

	if(read_commons(ini, com)==-1)
		return -1;

	if(read_channel_list(ini,det.boloname)==-1)
		return -1;

	if(read_fits_file_list(ini, dir,samples_struct)==-1)
		return -1;

	if(read_user_options(ini,u_opt)==-1)
		return -1;

	/*if(read_noise_file_list(ini, extentnoiseSP)==-1)
		return -1;*/

	if(read_ell_file(ini, ellFile)==-1)
		return -1;


	if(read_map_file(ini, signame)==-1)
		return -1;

	if(read_mixmatfile(ini, MixMatfile)==-1)
		return -1;


	printf("\nsanePS parser operations completed :\n");
	cout << "You have specified the following options : \n\n";

	print_directories(dir);
	print_commons(com);
	print_parser(u_opt);



	samples_struct.ntotscan = (samples_struct.fitsvect).size();
	det.ndet = det.boloname.size();	// ndet = number of detectors

	// Check improper usage
	if (det.ndet== 0) {
		cerr << "Must provide at least one channel.\n\n";
		return -1;
		//usage(argv[0]);
	}

	if(samples_struct.ntotscan == 0){
		cerr << "Must provide at least one scan.\n\n";
		return -1;
	}

	printf("Number of scans      : %ld\n",samples_struct.ntotscan);
	printf("Number of bolometers : %ld\n",det.ndet);

	//nnf = (int)extentnoiseSP.size();
	//nnf=1; // Temporarily
	/*if (nnf != 1 && nnf != samples_struct.ntotscan){
		cerr << "ERROR: There should be one noise power spectrum file per scan, or a single one for all the scans. Check -K options" << endl;
		exit(1);
	}

	if (nnf == 1 && samples_struct.ntotscan > 1)
		extentnoiseSP.resize(samples_struct.ntotscan, extentnoiseSP[0]);*/

	//  printf("%d\n",nnf);



	iniparser_freedict(ini);
	return 0 ;
}
