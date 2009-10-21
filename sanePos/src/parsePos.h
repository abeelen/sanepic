/*
 * parsePos.h
 *
 *  Created on: 11 juin 2009
 *      Author: matthieu
 */

#ifndef PARSEPOS_H_
#define PARSEPOS_H_

#include <vector>
#include <string>

#include "mpi_architecture_builder.h"

using namespace std;
/*!
 * - Parse sanePos input options using the ini file given in the command line \n
 * - Verify the good usage of each option, deal with incorrect values, warn the user when a line is \n
 * missing in the ini file !
 * - Initialise sanePos variable
 */
int parse_sanePos_ini_file(char * ini_name,struct user_options_sanepos &u_opt,
		long &ntotscan, long &ndet,
		std::vector<string> &bolonames, long *&nsamples,
		std::vector<struct box> &boxFile, std::vector<string> &fitsvect, std::vector<long> &scans_index);

#endif /* PARSEPOS_H_ */
