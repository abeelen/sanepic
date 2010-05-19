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
int parse_saneInv_ini_file(char * ini_name,struct samples &samples_struct,struct common &dir,std::string &boloname,  int rank);

#endif /* PARSEINV_H_ */
