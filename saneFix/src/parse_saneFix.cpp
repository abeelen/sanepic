

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
#include "tools.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#include "inputFileIO.h"
#include "parse_saneFix.h"

using namespace std;


int parse_saneFix_ini_file(char * ini_name, struct common &dir,
		struct samples &samples_struct, int rank)
{


	dictionary	*	ini ;


	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	// get directories path
	if(read_common(ini, dir, rank)==1)
		return -1;

	// get fits file that have to be fixed
	if(read_fits_file_list(ini, dir,samples_struct, rank)==1)
		return -1;

	// store number of scans
	samples_struct.ntotscan = (samples_struct.fitsvect).size();

	// free iniparser dictionnary
	iniparser_freedict(ini);
	return 0;
}
