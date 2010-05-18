

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


#include "inputFileIO.h"
#include "parse_saneMerge.h"

using namespace std;


int parse_saneMerge_ini_file(char * ini_name[], int arg, std::string &dir,
		struct samples &samples_struct)
/*! Parse user command line */
{

	dir=ini_name[1]; // ini file
	if(check_path(dir, "saneMerge Output Directory")) // get output path from command line, check its validity
		return 2;

	samples_struct.ntotscan = 0;
	for(int ii=2;ii<arg;ii++){ // get input file names from command line
		samples_struct.fitsvect.push_back(ini_name[ii]);
		samples_struct.ntotscan++;
	}

	readFrames(samples_struct.fitsvect, samples_struct.nsamples); // for each file, read and store number of samples in nsamples tab

	return 0;

}


