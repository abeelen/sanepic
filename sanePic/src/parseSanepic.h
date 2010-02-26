/*
 * parseSanepic.h
 *
 *  Created on: 18 juin 2009
 *      Author: matthieu
 */

#ifndef PARSESANEPIC_H_
#define PARSESANEPIC_H_


#include <vector>
#include "mpi_architecture_builder.h"


/*!
 * - Parse sanePic input options using the ini file given in the command line \n
 * - Verify the good usage of each option, deal with incorrect values, warn the user when a line is \n
 * missing in the ini file !
 * - Initialise sanePic variable
 */


int parse_sanePic_ini_file(char * ini_name,struct param_process &proc_param, struct param_positions &pos_param, int &iterw, struct common &dir, struct samples &samples_struct,
		struct detectors &det, std::vector<double> &fcut, int rank);

#endif /* PARSESANEPIC_H_ */
