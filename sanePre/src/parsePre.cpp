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
#include "covMatrix_IO.h"

extern "C"{
#include "nrutil.h"
#include "iniparser.h"
#include "dictionary.h"
}

using namespace std;

int parse_sanePre_ini_file(char * ini_name,struct param_process &proc_param, struct param_positions &pos_param, struct common &dir, struct samples &samples_struct,
		struct detectors &det, std::vector<double> &fcut,int rank, int size)

{
	dictionary	*	ini ;

	//int nnf=0; /*! extentnoiseSp_list number of elements*/


	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) { // if dictionnary was not found, return an error message and exit
		if(rank==0)
			fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return 2 ;
	}

	//DEFAULT PARAMETERS
	proc_param.napod = 0; /*! number of samples to apodize*/
	proc_param.fsamp = 0.0;// 25.0; /*! sampling frequency : BLAST Specific*/

	//Parser parameter (Program options)
	proc_param.NORMLIN = 0; /*!  baseline is removed from the data, NORMLIN = 1 else 0 */
	proc_param.NOFILLGAP = 0; /*! fill the gap ? default is YES*/
	proc_param.CORRon = 1; /*! correlation included in the analysis (=1), else 0, default 0*/
	proc_param.remove_polynomia = 1; /*! remove a polynomia fitted to the data*/
	proc_param.f_lp = 0.0; // low pass filter frequency
	pos_param.flgdupl = 0; // map duplication factor


	samples_struct.ntotscan=0; /*! total number of scans*/
	det.ndet=0; /*! number of channels*/



	// printf dictionnary to stderr for debugging
	//iniparser_dump(ini, stderr);


	/* Get sanepic_preprocess attributes */


	if(read_common(ini, dir,rank))
		return 2;

	if(read_channel_list(ini,dir,det.boloname,rank))
		return 2;

	if(read_fits_file_list(ini, dir,samples_struct,rank))
		return 2;

	if(read_param_process(ini,proc_param,rank))
		return 2;

	if(read_param_positions(ini,pos_param,rank))
		return 2;

	//if(read_noise_file_list(ini, extentnoiseSP)==-1)
	//return -1;

	if(read_noise_cut_freq(ini, proc_param, fcut,rank))
		return 2;

	if(rank==0){
		//		printf("\nsanePre parser operations completed :\n");
		cout << "You have specified the following options : \n\n";

		print_common(dir);
		print_param_process(proc_param);
		print_param_positions(pos_param);
	}

	samples_struct.ntotscan = (samples_struct.fitsvect).size();
	det.ndet = det.boloname.size();	// ndet = number of detectors

	if(size>det.ndet)
		return 3;

	// Check improper usage
	if (det.ndet== 0) {
		if(rank==0){
			cerr << "Must provide at least one channel.\n\n";
		}
		return 2;
	}

	if(samples_struct.ntotscan == 0){
		if(rank==0)
			cerr << "Must provide at least one scan.\n\n";
		return 2;
	}

	if(rank==0){
		printf("Number of scans      : %ld\n",samples_struct.ntotscan);
		printf("Number of bolometers : %ld\n",det.ndet);
	}


	// the number of noise cutting frequency must be egal to one (same for all scans) or ntotscan (one per scan)
	//	if ((int)fcut.size()==0){
	//		cerr << "Please give a correct number of noise cut frequency : 1 or 1 per scan\n";
	//		return 2;
	//	}

	if((int)fcut.size()==0){
		string fname;
		std::vector<string> bolos;
		long nbins;
		double *ell;
		double **Rellth;

		read_cov_matrix_file(ini, fname, rank);
		fname = dir.noise_dir + fname;
		read_CovMatrix(fname, bolos, nbins, ell, Rellth);
		fcut.push_back(ell[0]);
//		cout << ell[0] << " " << fcut[0] << endl;
		delete [] ell;
		long nBolos=bolos.size();
		free_dmatrix(Rellth, 0, nBolos * nBolos - 1, 0, nbins - 1);
	}

	// if only one fcut, extend to all scans
	if((int)fcut.size()==1)
		fcut.resize(samples_struct.ntotscan, fcut[0]);


	// cleaning up
	iniparser_freedict(ini);


	return 0 ;
}
