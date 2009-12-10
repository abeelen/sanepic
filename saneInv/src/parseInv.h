/*
 * parseInv.h
 *
 *  Created on: 18 juin 2009
 *      Author: matthieu
 */

#ifndef PARSEINV_H_
#define PARSEINV_H_


#include <string>

/*!
 * - Parse saneInv input options using the ini file given in the command line \n
 * - Verify the good usage of each option, deal with incorrect values, warn the user when a line is \n
 * missing in the ini file !
 * - Initialise saneInv variable
 */
int parse_saneInv_ini_file(char * ini_name, std::string &fname,struct samples &samples_struct,struct directories &dir,std::string &boloname, std::string &noiseSp_dir_output, std::string &extentnoiseSp);
std::string		Basename(std::string path);
//std::string		Basename(char * path);
#endif /* PARSEINV_H_ */
