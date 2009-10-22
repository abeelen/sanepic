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
int parse_sanePos_ini_file(char * ini_name,struct input_commons &com, struct directories &dir,
		struct detectors &det,struct samples &samples_struct,
		std::vector<struct box> &boxFile);

#endif /* PARSEPOS_H_ */
