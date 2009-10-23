/*
 * parseSanepic.h
 *
 *  Created on: 18 juin 2009
 *      Author: matthieu
 */

#ifndef PARSESANEPIC_H_
#define PARSESANEPIC_H_


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <list>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include "dataIO.h"


extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

using namespace std;

/*!
 * - Parse sanePic input options using the ini file given in the command line \n
 * - Verify the good usage of each option, deal with incorrect values, warn the user when a line is \n
 * missing in the ini file !
 * - Initialise sanePic variable
 */

/*int parse_sanePic_ini_file(char * ini_name,struct user_options &u_opt, int &iterw, long &ntotscan, long &ndet,
		std::vector<string> &bolonames,long *&nsamples,
		std::vector<struct box> & boxFile, std::vector<double> &fcut, std::vector<string> &extentnoiseSP,std::vector<string> &fitsvect,std::vector<string> &noisevect, std::vector<long> &scans_index);
*/

int parse_sanePic_ini_file(char * ini_name,struct user_options &u_opt, int &iterw, struct directories &dir, struct samples &samples_struct,struct input_commons &com,
		struct detectors &det, std::vector<struct box> & boxFile,std::vector<string> &extentnoiseSP, std::vector<double> &fcut);

#endif /* PARSESANEPIC_H_ */
