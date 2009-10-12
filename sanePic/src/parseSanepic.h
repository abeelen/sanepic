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
#include "positionsIO.h"


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
/*int parse_sanePic_ini_file(char * ini_name, double &pixdeg, int  &shift_data_to_point, long  &napod,double &fsamp, bool &NOFILLGAP,bool &NORMLIN,bool &projgaps,bool &remove_polynomia, bool &flgdupl,
		bool &CORRon, int &iterw, long &ntotscan, long &ndet, double &f_lp, double &f_lp_Nk, string &dirfile, string &outdir, string &tmp_dir, string &termin,
		string &MixMatfile, std::vector<string> &bolonames,long *&fframes,long *&nsamples, string &fname,
		std::vector<long> &xxi, std::vector<long> &xxf, std::vector<long> &yyi, std::vector<long> &yyf, std::vector<double> &fcut, std::vector<string> &extentnoiseSP,std::vector<string> &fitsvect,std::vector<string> &noisevect, std::vector<long> &scans_index);*/

int parse_sanePic_ini_file(char * ini_name,struct user_options &u_opt, int &iterw, long &ntotscan, long &ndet,
		string &MixMatfile, std::vector<string> &bolonames,long *&fframes,long *&nsamples,
		std::vector<struct box> & boxFile, std::vector<double> &fcut, std::vector<string> &extentnoiseSP,std::vector<string> &fitsvect,std::vector<string> &noisevect, std::vector<long> &scans_index);

#endif /* PARSESANEPIC_H_ */
