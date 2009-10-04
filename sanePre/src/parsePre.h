/*
 * parsePre.h
 *
 *  Created on: 12 juin 2009
 *      Author: matthieu
 */

#ifndef PARSEPRE_H_
#define PARSEPRE_H_

#include <vector>
#include <string>

using namespace std;

/*!
 * - Parse sanePre input options using the ini file given in the command line \n
 * - Verify the good usage of each option, deal with incorrect values, warn the user when a line is \n
 * missing in the ini file !
 * - Initialise sanePre variable
 */

/*int parse_sanePre_ini_file(char * ini_name, int  &shift_data_to_point, long  &napod,double &fsamp, bool &NOFILLGAP,bool &NORMLIN,bool &remove_polynomia, bool &flgdupl,
		bool &CORRon, long &ntotscan, long &ndet, double &f_lp, string &dirfile, string &outdir,
		 string &noiseSppreffile,
		int &coordsyst, std::vector<string> &bolonames,long *&fframes, long *&nsamples,std::vector<struct box> & boxFile,
		std::vector<string> &extentnoiseSP, std::vector<double> &fcut,std::vector<string> &fitsvect,std::vector<string> &noisevect, std::vector<long> &scans_index);*/

int parse_sanePre_ini_file(char * ini_name,struct user_options &u_opt,
		long &ntotscan, long &ndet,
		std::vector<string> &bolonames,long *&fframes, long *&nsamples,
		std::vector<struct box> & boxFile,std::vector<string> &extentnoiseSP, std::vector<double> &fcut, std::vector<string> &fitsvect,std::vector<string> &noisevect, std::vector<long> &scans_index);


#endif /* PARSEPRE_H_ */
