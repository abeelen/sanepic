/*
 * parse_saneCheck.h
 *
 *  Created on: 3 déc. 2009
 *      Author: matthieu
 */

#ifndef PARSE_SANECHECK_H_

int parse_saneCheck_ini_file(char * ini_name, struct common &dir,
		std::vector<detectors> detector_tab,struct samples &samples_struct, double &fsamp, int rank);


#define PARSE_SANECHECK_H_


#endif /* PARSE_SANECHECK_H_ */
