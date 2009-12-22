/*
 * parseSanepic.cpp
 *
 *  Created on: 18 juin 2009
 *      Author: matthieu
 */

#include <iostream>

#include "parseSanepic.h"
#include "parser_functions.h"
#include "struct_definition.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

using namespace std;

int parse_sanePic_ini_file(char * ini_name,struct param_process &proc_param, struct param_positions &pos_param, int &iterw, struct directories &dir, struct samples &samples_struct,
		struct detectors &det, std::vector<double> &fcut, int rank)
{

	dictionary	*	ini ;

	//int nnf=0; /*! extentnoiseSp_list number of elements*/

	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}



	//DEFAULT PARAMETERS
	proc_param.napod = 0; /*!  number of samples to apodize */
	proc_param.fsamp = 0.0; //25.0; /*!  sampling frequency : BLAST Specific */
	iterw = 10; /*!  period in iterations to which the data are written to disk, 0 = no intermediate map to be written*/
	proc_param.NORMLIN = 0; /*!  baseline is removed from the data, NORMLIN = 1 else 0 */
	proc_param.NOFILLGAP = 0; /*!  fill the gap ? default is YES (debug parameter) */
	proc_param.remove_polynomia = 1; /*! Remove a fitted polynomia from the data ? */
	proc_param.CORRon = 1; /*!  correlation included in the analysis (=1), else 0, default 0 */
	proc_param.f_lp=0.0; /*! frequencies : filter knee freq, noise PS threshold freq ; frequencies converted in a number of samples */
	pos_param.flgdupl = 0; // map duplication boolean

	samples_struct.ntotscan=0; /*! total number of scans */
	det.ndet=0; /*! number of channels */

	// printf dictionnary to stderr for debugging
	//iniparser_dump(ini, stderr);

	printf("\nsanepic_conjugate_gradient:\n");

	if(read_directories(ini, dir, rank)==-1)
		return -1;

	if(read_param_process(ini, proc_param, rank)==-1)
		return -1;

	if(read_param_positions(ini, pos_param, rank)==-1)
		return -1;

	if(read_channel_list(ini,det.boloname, rank)==-1)
		return -1;

	if(read_fits_file_list(ini, dir,samples_struct, rank)==-1)
		return -1;



	/*if(read_noise_file_list(ini, extentnoiseSP)==-1)
		return -1;*/

	if(read_noise_cut_freq(ini, fcut, rank)==-1)
		return -1;

	if(read_iter(ini, iterw, rank)==-1)
		return -1;

	if(rank==0){

		printf("\nsanePre parser operations completed :\n");
		cout << "You have specified the following options : \n\n";

		print_directories(dir);
		print_param_process(proc_param);
		print_param_process(proc_param);
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


	// the number of noise cutting frequency must be egal to one (same for all scans) or ntotscan (one per scan)
	if ((int)fcut.size()==0){
		cerr << "Please give a correct number of noise cut frequency : 1 or 1 per scan\n";
		exit(0);
	}

	// if only one fcut, extend to all scans
	if(fcut.size()==1)
		fcut.resize(samples_struct.ntotscan, fcut[0]);


	// cleaning up
	iniparser_freedict(ini);

	return 0 ;
}
