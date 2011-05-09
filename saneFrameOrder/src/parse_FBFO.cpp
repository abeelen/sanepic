#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <string>
#include <unistd.h>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>

#include "dataIO.h"
#include "mpi_architecture_builder.h"
#include "parser_functions.h"
#include "struct_definition.h"
#include "inputFileIO.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#include "parse_FBFO.h"

using namespace std;

int parse_FBFO(char * ini_name, string &parser_output, struct samples &samples_struct,struct param_common &dir)
{
	dictionary	*	ini ;

	string str, dirfile;

	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	samples_struct.ntotscan=0; /*! total number of scans*/

	read_common(parser_output, ini, dir);

	// Fill fitsvec, noisevect, scans_index with values read from the 'str' filename
	if(read_fits_list(parser_output, dir.input_dir+dir.fits_filelist, samples_struct)!=0)
		return -1;

	cout << "You have specified the following options : \n\n";

	print_common(dir);

	samples_struct.ntotscan = (samples_struct.fitsvect).size();

	if(samples_struct.ntotscan == 0){
		cerr << "Must provide at least one scan.\n\n";
		return -1;
	}

	printf("Number of scans      : %ld\n",samples_struct.ntotscan);

	iniparser_freedict(ini);

	return 0 ;
}


