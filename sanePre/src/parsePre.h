/*
 * parsePre.h
 *
 *  Created on: 12 juin 2009
 *      Author: matthieu
 */

#ifndef PARSEPRE_H_
#define PARSEPRE_H_

#include <vector>
#include "mpi_architecture_builder.h"



/*!
 * - Parse sanePre input options using the ini file given in the command line \n
 * - Verify the good usage of each option, deal with incorrect values, warn the user when a line is \n
 * missing in the ini file !
 * - Initialise sanePre variable
 */

int parse_sanePre_ini_file(char * ini_name,struct param_process &proc_param,struct param_positions &pos_param, struct common &dir, struct samples &samples_struct,
		struct detectors &det, std::vector<double> &fcut, int rank, int size);

#endif /* PARSEPRE_H_ */
