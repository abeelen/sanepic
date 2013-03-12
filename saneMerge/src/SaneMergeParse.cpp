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

#include "DataIO.h"
#include "StructDefinition.h"
#include "MPIConfiguration.h"
#include "ParserFunctions.h"
#include "SaneMergeTools.h"
#include "InputFileIO.h"
#include "SaneMergeParse.h"

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

	if(check_common(output, dir, 0))
		return 2;

	samples_struct.ntotscan = 0;
	for(int ii=2;ii<arg;ii++){ // get input file names from command line
		samples_struct.fitsvect.push_back(opt_name[ii]);
		samples_struct.ntotscan++;
	}


	cout << "You have specified the following options : \n\n";

	print_common(dir); /* print dir locations on stdout */

	iniparser_freedict(ini);

	return 0;

}


