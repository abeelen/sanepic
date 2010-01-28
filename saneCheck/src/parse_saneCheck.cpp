

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

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#include "inputFileIO.h"
#include "parse_saneCheck.h"

using namespace std;


int parse_saneCheck_ini_file(char * ini_name, struct directories &dir,
		struct detectors &det,struct samples &samples_struct, int rank)
{


	dictionary	*	ini ;

	struct param_positions pos_param;
	struct param_process proc_param;
	std::vector<double> fcut;
	double fcut_double;
	string MixMatfile, ellFile, signame;
	long ncomp;

	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	if(read_directories(ini, dir, rank)==1)
		return -1;

	if(read_channel_list(ini,det.boloname, rank)==1)
		return -1;

	if(read_fits_file_list(ini, dir,samples_struct, rank)==1)
		return -1;

	samples_struct.ntotscan = (samples_struct.fitsvect).size();
	det.ndet = (long)((det.boloname).size());

	if (det.ndet == 0) {
		cerr << "Must provide at least one channel.\n\n";
		return -1 ;
	}

	if(read_param_positions(ini, pos_param, rank)==1)
		return -1;

	if(read_param_process(ini, proc_param, rank)==1)
		return -1;

	if(read_noise_cut_freq(ini, fcut,rank)==1)
		return -1;

	if(read_ell_file(ini, ellFile, rank) ||
			read_map_file(ini, signame, rank) ||
			read_mixmatfile(ini, MixMatfile, rank)||
			read_ncomp(ini, ncomp, rank) ||
			read_fcut(ini, fcut_double, rank))
		return -1;



	if(rank==0){

		//		printf("\nsaneCheck parser operations completed :\n");
		cout << "You have specified the following options : \n\n";

		print_directories(dir);
		print_param_process(proc_param);
		print_param_positions(pos_param);


		printf("Number of scans      : %ld\n",samples_struct.ntotscan);
		printf("Number of bolometers : %ld\n",det.ndet);
	}


	iniparser_freedict(ini);
	return 0 ;
}
