/*
 * parse_FBFO.h
 *
 *  Created on: 15 juin 2009
 *      Author: matthieu
 */

#ifndef PARSE_FBFO_H_
#define PARSE_FBFO_H_

//#include <string>

//int parse_FBFO(char * ini_name, std::string &tmp_dir, long &ntotscan, long *&nsamples,
//		std::vector<std::string> &fitsvect,std::vector<std::string> &noisevect,std::vector<int> &scans_index, int rank);
int parse_FBFO(char * ini_name,struct samples &samples_struct,struct common &dir);

#endif /* PARSE_FBFO_H_ */
