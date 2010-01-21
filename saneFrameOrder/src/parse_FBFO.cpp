/*
 * parse_FBFO.cpp
 *
 *  Created on: 15 juin 2009
 *      Author: matthieu
 */

/*
 * testPos.cpp
 *
 *  Created on: 11 juin 2009
 *      Author: matthieu
 */

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

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#include "parse_FBFO.h"

using namespace std;

//int parse_FBFO(char * ini_name, string &tmp_dir, long &ntotscan, long *&nsamples,
//		std::vector<string> &fitsvect, std::vector<string> &noisevect, std::vector<int> &scans_index, int rank)
int parse_FBFO(char * ini_name,struct samples &samples_struct,struct directories &dir)
{
	dictionary	*	ini ;

	/* Some temporary variables to hold query results */
//	char		*	s ;
	string str, dirfile;

	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	// printf dictionnary to stderr for debugging
	//iniparser_dump(ini, stderr);


	samples_struct.ntotscan=0; /*! total number of scans*/

	if(read_directories(ini, dir,0)==-1)
		return -1;

	if(read_fits_file_list(ini, dir,samples_struct,0)==-1)
		return -1;


	printf("\nsaneFrameOrder parser operations completed :\n");
	cout << "You have specified the following options : \n\n";

	print_directories(dir);



	samples_struct.ntotscan = (samples_struct.fitsvect).size();

	if(samples_struct.ntotscan == 0){
		cerr << "Must provide at least one scan.\n\n";
		return -1;
	}

	printf("Number of scans      : %ld\n",samples_struct.ntotscan);

	iniparser_freedict(ini);

	return 0 ;
}


