
#include <iostream>
#include <sstream>
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
#include "struct_definition.h"
#include "mpi_architecture_builder.h"
#include "parser_functions.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#include "inputFileIO.h"
#include "parse_saneSplit.h"

using namespace std;

int parse_saneSplit_ini_file(char * ini_name, struct common &dir)
{


	dictionary	*	ini ;


	// load dictionnary
	ini = iniparser_load(ini_name);


	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return 2 ;
	}

	if(read_common(ini, dir, 0)==1)
		return 2;



	//	if(read_channel_list(ini,det.boloname, 0)==1)
	//		return 2;

	//if(read_fits_file_list(ini, dir,samples_struct, 0)==1)
	//return 2;

	//	samples_struct.ntotscan = (samples_struct.fitsvect).size();
	//	det.ndet = (long)((det.boloname).size());
	//
	//	if (det.ndet == 0) {
	//		cerr << "Must provide at least one channel.\n\n";
	//		return 2;
	//	}

	//	if(read_parser_string(ini, "saneSplit:cut_limits_file", 0, fname)==1)
	//		return 2;



	cout << "You have specified the following options : \n\n";

	print_common(dir);
	//	printf("Number of scans      : %ld\n",samples_struct.ntotscan);
	//	printf("Number of bolometers : %ld\n",det.ndet);

	return 0;
}





