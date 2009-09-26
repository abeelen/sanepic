/*
 * parsePre.h
 *
 *  Created on: 12 juin 2009
 *      Author: matthieu
 */

#ifndef PARSEPRE_H_
#define PARSEPRE_H_

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
//int parse_sanePos_ini_file(char * ini_name);
int parse_sanePre_ini_file(char * ini_name, int  &shift_data_to_point, long  &napod,double &fsamp, bool &NOFILLGAP,bool &NORMLIN, bool &flgdupl,
		bool &CORRon, long &ntotscan, long &ndet, int &nnf,	double &f_lp, double &f_lp_Nk, string &dirfile, string &outdir, string &poutdir, string &bextension,
		string &fextension, string &pextension, string &termin, string &noiseSppreffile,
		int &coordsyst, std::vector<string> &bolonames,std::vector<long> &fframes_vec, std::vector<long> &nsamples_vec, string &fname, double &pixdeg);


#endif /* PARSEPRE_H_ */
