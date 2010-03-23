/*
 * parsePos.cpp
 *
 *  Created on: 11 juin 2009
 *      Author: matthieu
 */


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <unistd.h>
#include <string>
#include <vector>
#include <algorithm>

#include "dataIO.h"
#include "mpi_architecture_builder.h"
#include "parser_functions.h"
#include "struct_definition.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#include "inputFileIO.h"
#include "parsePos.h"

using namespace std;

/*int parse_sanePos_ini_file(char * ini_name,struct param_process_sanepos &proc_param,
		long &ntotscan, long &ndet,
		std::vector<string> &bolonames, long *&nsamples,
		std::vector<struct box> &boxFile, std::vector<string> &fitsvect, std::vector<long> &scans_index)*/
int parse_sanePos_ini_file(char * ini_name,struct param_process &proc_param, struct param_positions &pos_param, struct common &dir,
		struct detectors &det,struct samples &samples_struct,
		int rank)
{

	dictionary	*	ini ;

	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "\ncannot parse file: %s\n", ini_name);
		return -1 ;
	}



	// default values :
	proc_param.napod  = 0; /*! number of samples to apodize, =0 -> no apodisation */
	proc_param.NOFILLGAP = 0; /*! dont fill the gaps ? default is NO => the program fill */
	samples_struct.ntotscan=0; /*! total number of scans */
	det.ndet=0; /*! number of channels used*/
	pos_param.flgdupl = 0; // map duplication factor

	pos_param.maskfile = "";

	// printf dictionnary to stderr for debugging
	//iniparser_dump(ini, stderr);


	/* Get sanepic_compute_positions attributes */

	if(read_param_positions(ini, pos_param, rank))
		return -1;

	if(read_common(ini, dir, rank))
		return -1;

	if(read_param_process(ini, proc_param, rank))
		return -1;

	if(read_channel_list(ini, dir, det.boloname, rank))
		return -1;

	if(read_fits_file_list(ini, dir, samples_struct, rank))
		return -1;

	if(check_dirfile_paths(dir.tmp_dir)){
		cout << "Exiting...\n";
		return -1;
	}


	if(rank==0){

		//		printf("\nsanePos parser operations completed :\n");
		cout << "\nYou have specified the following options : \n";

		print_common(dir);
		print_param_process(proc_param);
		print_param_positions(pos_param);
	}

	/*
	if (tmpcount == 1 || tmpcount == 2 || tmpcount == 3){
		cerr << "ERROR: None or all the following keywords must be set: RA_min RA_max DEC_min DEC_max. Exiting. \n";
		return -1 ;
	}
	if (tmpcount2 == 1 || tmpcount2 == 2){
		cerr << "ERROR: None or all the following keywords must be set: RA_source DEC_source map_radius. Exiting. \n";
		return -1 ;
	}



	//if (tmpcount == 4)
		//proc_param.bfixc = 1;

	if (tmpcount2 == 3){
		if (tmpcount == 4){
			cerr << "ERROR: Conflicting input parameter: RA_min RA_max DEC_min DEC_max keywords are not compatible with RA_source DEc_source map_radius keywords . Exiting. \n";
			return -1 ;
		}
		//proc_param.bfixc = 1;
		proc_param.coordscorner[0] = proc_param.srccoord[0];
		proc_param.coordscorner[1] = proc_param.srccoord[0];
		proc_param.coordscorner[2] = proc_param.srccoord[1];
		proc_param.coordscorner[3] = proc_param.srccoord[1];
	}

	 */

	samples_struct.ntotscan = (samples_struct.fitsvect).size();
	det.ndet = (long)((det.boloname).size());

	if (det.ndet == 0) {
		cerr << "Must provide at least one channel.\n\n";
		return -1 ;
		//usage(argv[0]);
	}

	//cout << "framegiven : " << samples_vct.framegiven << endl;

	if(rank==0){
		printf("Number of scans      : %ld\n",samples_struct.ntotscan);
		printf("Number of bolometers : %ld\n",det.ndet);
	}

	iniparser_freedict(ini);
	return 0 ;
}



