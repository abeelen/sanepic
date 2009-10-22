/*
 * parsePos.cpp
 *
 *  Created on: 11 juin 2009
 *      Author: matthieu
 */


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

#include "positionsIO.h"
#include "mpi_architecture_builder.h"
#include "parser_functions.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#include "inputFileIO.h"
#include "parsePos.h"

/*int parse_sanePos_ini_file(char * ini_name,struct user_options_sanepos &u_opt,
		long &ntotscan, long &ndet,
		std::vector<string> &bolonames, long *&nsamples,
		std::vector<struct box> &boxFile, std::vector<string> &fitsvect, std::vector<long> &scans_index)*/
int parse_sanePos_ini_file(char * ini_name,struct input_commons &com, struct directories &dir,
		long &ndet,	std::vector<string> &bolonames,struct samples &samples_str,
		std::vector<struct box> &boxFile, struct samples_vect &samples_vct)
{

	dictionary	*	ini ;

	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}


	// printf dictionnary to stderr for debugging
	//iniparser_dump(ini, stderr);


	/* Get sanepic_compute_positions attributes */
	printf("sanepic_compute_positions:\n");



	if(read_dirfile(ini, dir)==-1)
		return -1;

	if(read_tmpdir(ini, dir)==-1)
		return -1;

	if(read_outdir(ini, dir)==-1)
		return -1;

	if(read_commons(ini, com)==-1)
		return -1;

	if(read_channel_list(ini,bolonames)==-1)
		return -1;

	if(read_fits_file_list(ini, dir,samples_vct,samples_str)==-1)
		return -1;

	if(read_box_coord(ini,boxFile)==-1)
		return -1;





	printf("\nsanePos parser operations completed :\n");
	cout << "You have specified the following options : \n";

	print_directories(dir);
	print_commons(com);

	/*
	if (tmpcount == 1 || tmpcount == 2 || tmpcount == 3){
		cerr << "ERROR: None or all the following keywords must be set: RA_min RA_max DEC_min DEC_max. Exiting. \n";
		return -1 ;
	}
	if (tmpcount2 == 1 || tmpcount2 == 2){
		cerr << "ERROR: None or all the following keywords must be set: RA_source DEC_source map_radius. Exiting. \n";
		return -1 ;
	}



	//if (tmpcount == 4)
		//u_opt.bfixc = 1;

	if (tmpcount2 == 3){
		if (tmpcount == 4){
			cerr << "ERROR: Conflicting input parameter: RA_min RA_max DEC_min DEC_max keywords are not compatible with RA_source DEc_source map_radius keywords . Exiting. \n";
			return -1 ;
		}
		//u_opt.bfixc = 1;
		u_opt.coordscorner[0] = u_opt.srccoord[0];
		u_opt.coordscorner[1] = u_opt.srccoord[0];
		u_opt.coordscorner[2] = u_opt.srccoord[1];
		u_opt.coordscorner[3] = u_opt.srccoord[1];
	}

*/

	samples_str.ntotscan = (samples_vct.fitsvect).size();
	ndet = (long)bolonames.size();

	if (ndet == 0) {
		cerr << "Must provide at least one channel.\n\n";
		return -1 ;
		//usage(argv[0]);
	}

	//cout << "framegiven : " << samples_vct.framegiven << endl;


	printf("Number of scans      : %ld\n",samples_str.ntotscan);
	printf("Number of bolometers : %ld\n",ndet);

	iniparser_freedict(ini);
	return 0 ;
}



