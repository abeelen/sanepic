/*
 * parsePre.cpp
 *
 *  Created on: 12 juin 2009
 *      Author: matthieu
 */


#include <iostream>


#include "parsePre.h"
#include "parser_functions.h"
#include "struct_definition.h"


extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

using namespace std;

/*
int parse_sanePre_ini_file(char * ini_name,struct user_options &u_opt, struct directories &dir, struct samples &samples_struct,struct input_commons &com,
		struct detectors &det, std::vector<struct box> & boxFile,std::vector<string> &extentnoiseSP, std::vector<double> &fcut)
 */
int parse_sanePre_ini_file(char * ini_name,struct user_options &u_opt, struct directories &dir, struct samples &samples_struct,struct input_commons &com,
		struct detectors &det, std::vector<struct box> & boxFile, std::vector<double> &fcut,int rank)

{
	dictionary	*	ini ;

	//int nnf=0; /*! extentnoiseSp_list number of elements*/


	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) { // if dictionnary was not found, return an error message and exit
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	// printf dictionnary to stderr for debugging
	//iniparser_dump(ini, stderr);


	/* Get sanepic_preprocess attributes */
	printf("\n[%d] sanepic_preprocess\n",rank);

	if(read_directories(ini, dir,rank)==-1)
		return -1;

	if(read_commons(ini, com,rank)==-1)
		return -1;

	if(read_channel_list(ini,det.boloname,rank)==-1)
		return -1;

	if(read_fits_file_list(ini, dir,samples_struct,rank)==-1)
		return -1;

	if(read_box_coord(ini,boxFile,rank)==-1)
		return -1;

	if(read_user_options(ini,u_opt,rank)==-1)
		return -1;

	//if(read_noise_file_list(ini, extentnoiseSP)==-1)
	//return -1;

	if(read_noise_cut_freq(ini, fcut,rank)==-1)
		return -1;

	if(rank==0){
		printf("\nsanePre parser operations completed :\n");
		cout << "You have specified the following options : \n\n";

		print_directories(dir);
		print_commons(com);
		print_parser(u_opt);
	}


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

	if(rank==0){
		printf("Number of scans      : %ld\n",samples_struct.ntotscan);
		printf("Number of bolometers : %ld\n",det.ndet);
	}



	//nnf = number of noise PS files
	//	nnf = (int)extentnoiseSP.size();

	//nnf=1; // Debug
	/*	if (nnf != 1 && nnf != samples_struct.ntotscan){
		cerr << "ERROR: There should be one noise power spectrum file per scan, or a single one for all the scans. Check -K options" << endl;
		return -1;
	}*/

	//if only one extension for all the noisePS file : extend to all the scans
	//if (nnf == 1 && samples_struct.ntotscan > 1)
	//extentnoiseSP.resize(samples_struct.ntotscan, extentnoiseSP[0]);

	// the number of noise cutting frequency must be egal to one (same for all scans) or ntotscan (one per scan)
	if ((int)fcut.size()==0){
		cerr << "Please give a correct number of noise cut frequency : 1 or 1 per scan\n";
		exit(0);
	}

	// if only one fcut, extend to all scans
	if((int)fcut.size()==1)
		fcut.resize(samples_struct.ntotscan, fcut[0]);


	// cleaning up
	iniparser_freedict(ini);


	return 0 ;
}
