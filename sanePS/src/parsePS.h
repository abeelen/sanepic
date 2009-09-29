/*
 * parsePS.h
 *
 *  Created on: 18 juin 2009
 *      Author: matthieu
 */


#ifndef PARSEPS_H_
#define PARSEPPS_H_

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
 * - Parse sanePS input options using the ini file given in the command line \n
 * - Verify the good usage of each option, deal with incorrect values, warn the user when a line is \n
 * missing in the ini file !
 * - Initialise sanePS variable
 */
int parse_sanePS_ini_file(char * ini_name, int  &shift_data_to_point, long  &napod,double &fsamp, bool &NOFILLGAP,bool &NORMLIN,bool &remove_polynomia, bool &flgdupl,
		long &ntotscan, long &ndet, string &dirfile, string &outdir, string &tmp_dir,/* string &bextension,
		string &fextension,*/ string &noiseSppreffile,
		std::vector<string> &bolonames,long *&fframes, long *&nsamples, std::vector<string> &extentnoiseSP, string &MixMatfile,string &signame,
		std::vector<string> &fitsvect,std::vector<string> &noisevect, std::vector<long> &scans_index);



#endif /* PARSEPPS_H_ */
