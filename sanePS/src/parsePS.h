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
/*int parse_sanePS_ini_file(char * ini_name, struct user_options &u_opt,
		long &ntotscan, long &ndet,
		std::vector<string> &bolonames, long *&nsamples, std::vector<string> &extentnoiseSP, string &MixMatfile, string & ellFile, string &signame,
		std::vector<string> &fitsvect,std::vector<string> &noisevect, std::vector<long> &scans_index);*/

int parse_sanePS_ini_file(char * ini_name, struct user_options &u_opt, struct directories &dir, struct samples &samples_struct,struct input_commons &com,
		struct detectors &det, std::string &MixMatfile, std::string &ellFile, std::string &signame);


#endif /* PARSEPPS_H_ */
