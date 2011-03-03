

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


int parse_saneFix_ini_file(char * ini_name, string &output, struct param_common &dir,
		struct samples &samples_struct, int rank)
{


	dictionary	*	ini ;


	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	default_param_common(dir);
	read_common(output, ini, dir);

	string filename;
	filename = dir.input_dir+dir.fits_filelist;
	if(read_fits_list(output, filename, samples_struct)!=0)
		return 1;
	samples_struct.ntotscan = (samples_struct.fitsvect).size();

	// No fundamental reason but we need to take fitsbasename anyway ! So we can add path here !
	for(int iframe=0;iframe<(int)((samples_struct.fitsvect).size());iframe++){
		samples_struct.fitsvect[iframe] = dir.dirfile + samples_struct.fitsvect[iframe];
	}

	readFrames(samples_struct.fitsvect, samples_struct.nsamples);


	// free iniparser dictionnary
	iniparser_freedict(ini);
	return 0;
}
