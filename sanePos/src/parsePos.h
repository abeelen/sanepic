/*
 * parsePos.h
 *
 *  Created on: 11 juin 2009
 *      Author: matthieu
 */

#ifndef PARSEPOS_H_
#define PARSEPOS_H_

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
#include "positionsIO.h"


extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

using namespace std;
/*!
 * - Parse sanePos input options using the ini file given in the command line \n
 * - Verify the good usage of each option, deal with incorrect values, warn the user when a line is \n
 * missing in the ini file !
 * - Initialise sanePos variable
 */
int parse_sanePos_ini_file(char * ini_name,bool &bfixc, int  &shift_data_to_point, long  &napod, bool &NOFILLGAP, bool &flgdupl,
		double * srccoord, double * coordscorner, double &radius, long &ntotscan, long &ndet, int &nnf,
		double &pixdeg, string &dirfile, string &outdir, string &poutdir, string &bextension,
		string &fextension, string &pextension, string &file_offsets, string &file_frame_offsets, /*string &termin,*/
		int &coordsyst, std::vector<string> &bolonames, std::vector<long> &fframes_vec, std::vector<long> &nsamples_vec,
		std::vector<long> &xxi,std::vector<long> &xxf, std::vector<long> &yyi, std::vector<long> &yyf);

#endif /* PARSEPOS_H_ */
