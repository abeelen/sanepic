

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

	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	if(read_directories(ini, dir, rank)==-1)
		return -1;

	if(read_channel_list(ini,det.boloname, rank)==-1)
		return -1;

	if(read_fits_file_list(ini, dir,samples_struct, rank)==-1)
		return -1;

	samples_struct.ntotscan = (samples_struct.fitsvect).size();
	det.ndet = (long)((det.boloname).size());

	if (det.ndet == 0) {
		cerr << "Must provide at least one channel.\n\n";
		return -1 ;
		//usage(argv[0]);
	}

	//cout << "framegiven : " << samples_vct.framegiven << endl;



	if(rank==0){

		printf("\nsaneCheck parser operations completed :\n");
		cout << "You have specified the following options : \n\n";

		print_directories(dir);


		printf("Number of scans      : %ld\n",samples_struct.ntotscan);
		printf("Number of bolometers : %ld\n",det.ndet);
	}


	iniparser_freedict(ini);
	return 0 ;
}
