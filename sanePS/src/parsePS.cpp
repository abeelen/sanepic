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
#include "struct_definition.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

using namespace std;

int parse_sanePS_ini_file(char * ini_name, struct param_process &proc_param, struct common &dir, struct samples &samples_struct,
		struct detectors &det,string &MixMatfile, string &ellFile, string &signame, int rank, long &ncomp, double &fcut)

{
	dictionary	*	ini ;

	//int nnf=0; // extentnoiseSp_list number of elements


	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		if(rank==0)
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}


	//DEFAULT PARAMETERS
	proc_param.napod = 0; // number of samples to apodize
	proc_param.fsamp = 0.0; //25.0; // sampling frequency : BLAST Specific
	proc_param.NORMLIN = 0; // baseline is removed from the data, NORMLIN = 1 else 0
	proc_param.NOFILLGAP = 0; // fill the gap ? default is YES (debug parameter)
	proc_param.remove_polynomia=1; // remove a fitted polynomia to the data ? 0 for no, > 0 = order of the poly

	samples_struct.ntotscan=0; // total number of scans
	det.ndet=0; // number of channels


	if( read_common(ini, dir, rank) ||
			read_param_process(ini, proc_param, rank) ||
			read_channel_list(ini,dir,det.boloname, rank) ||
			read_fits_file_list(ini, dir,samples_struct, rank) ||
			read_ell_file(ini, ellFile, rank) ||
			read_map_file(ini, signame, rank) ||
			read_mixmatfile(ini, MixMatfile, rank)||
			read_ncomp(ini, ncomp, rank) ||
			read_fcut(ini, fcut, rank) )
		exit(1);

	if(rank==0){
//		printf("\nsanePS parser operations completed :\n");
		cout << "\nYou have specified the following options : \n\n";

		print_common(dir);
		print_param_process(proc_param);

	}


	samples_struct.ntotscan = (samples_struct.fitsvect).size();
	det.ndet = det.boloname.size();	// ndet = number of detectors

	// Check improper usage
if (det.ndet < 2) {
		cerr << "Must provide at least two channels.\n\n";
		return -1;
		//usage(argv[0]);
	}

	if(samples_struct.ntotscan == 0){
		if(rank==0)
		cerr << "Must provide at least one scan.\n\n";
		return -1;
	}

	if(rank==0){
		printf("Number of scans      : %ld\n",samples_struct.ntotscan);
		printf("Number of bolometers : %ld\n",det.ndet);
	}


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
