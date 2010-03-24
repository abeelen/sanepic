

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
{

	dir=ini_name[1];
	if(check_path(dir, "saneMerge Output Directory"))
		return 2;

	samples_struct.ntotscan = 0;
	for(int ii=2;ii<arg;ii++){
		samples_struct.fitsvect.push_back(ini_name[ii]);
		samples_struct.ntotscan++;
	}

	readFrames(samples_struct.fitsvect, samples_struct.nsamples);

	return 0;

}


