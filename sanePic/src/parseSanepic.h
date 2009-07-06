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
#include "boloIO.h"


extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

using namespace std;
//int parse_sanePos_ini_file(char * ini_name);
int parse_sanePic_ini_file(char * ini_name, double &pixdeg, int  &shift_data_to_point, long  &napod,double &fsamp, bool &NOFILLGAP,bool &NORMLIN,bool &projgaps,bool &remove_polynomia, bool &flgdupl,
		bool &CORRon, int &iterw, bool &doInitPS, long &ntotscan, long &ndet, int &nnf,	double &f_lp, double &f_lp_Nk, string &dirfile, string &outdir, string &poutdir, string &bextension,
		string &fextension, string &pextension, string &termin, string &noiseSppreffile,
		int &coordsyst, string &MixMatfile, std::vector<string> &bolonames,std::vector<long> &fframes_vec, std::vector<long> &nsamples_vec, string &fname,
		std::vector<long> &xxi, std::vector<long> &xxf, std::vector<long> &yyi, std::vector<long> &yyf, std::vector<double> &fcut, std::vector<string> &extentnoiseSP);



#endif /* PARSESANEPIC_H_ */