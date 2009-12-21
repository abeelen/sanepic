/*
 * parsePS.h
 *
 *  Created on: 18 juin 2009
 *      Author: matthieu
 */


#ifndef PARSEPS_H_
#define PARSEPPS_H_



#include <string>
#include "mpi_architecture_builder.h"

/*!
 * - Parse sanePS input options using the ini file given in the command line \n
 * - Verify the good usage of each option, deal with incorrect values, warn the user when a line is \n
 * missing in the ini file !
 * - Initialise sanePS variable
 */
int parse_sanePS_ini_file(char * ini_name, struct param_process &proc_param, struct directories &dir, struct samples &samples_struct,
		struct detectors &det, std::string &MixMatfile, std::string &ellFile, std::string &signame, int rank, long &ncomp, double &fcut);


#endif /* PARSEPPS_H_ */
