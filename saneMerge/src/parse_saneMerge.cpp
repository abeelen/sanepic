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


int parse_saneMerge_ini_file(char * opt_name[], string &output, int arg, struct param_common &dir,
		struct samples &samples_struct)
/* Parse user command line */
{

	dictionary	*	ini ;


	// load dictionnary
	ini = iniparser_load(opt_name[1]);


	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", opt_name[1]);
		return 2 ;
	}

	// get directories path
	read_common(output, ini, dir);

	if(check_common(output, dir))
		return 2;

	samples_struct.ntotscan = 0;
	for(int ii=2;ii<arg;ii++){ // get input file names from command line
		samples_struct.fitsvect.push_back(opt_name[ii]);
		samples_struct.ntotscan++;
	}

	readFrames(dir.data_dir, samples_struct.fitsvect, samples_struct.nsamples); // for each file, read and store number of samples in nsamples tab

	cout << "You have specified the following options : \n\n";

	print_common(dir); /* print dir locations on stdout */

	iniparser_freedict(ini);

	return 0;

}


