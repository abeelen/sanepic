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

/*int parse_sanePre_ini_file(char * ini_name,struct user_options &u_opt,
		long &ntotscan, long &ndet,
		std::vector<string> &bolonames, long *&nsamples,
		std::vector<struct box> & boxFile,std::vector<string> &extentnoiseSP, std::vector<double> &fcut, std::vector<string> &fitsvect,std::vector<string> &noisevect, std::vector<long> &scans_index);
*/

int parse_sanePre_ini_file(char * ini_name,struct user_options &u_opt, struct directories &dir, struct samples &samples_struct,struct input_commons &com,
		struct detectors &det, std::vector<struct box> & boxFile, std::vector<double> &fcut);

#endif /* PARSEPRE_H_ */
