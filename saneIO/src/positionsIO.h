/*
 * positionsIO.h
 *
 *  Created on: 31 août 2009
 *      Author: abeelen
 */

#ifndef POSITIONSIO_H_
#define POSITIONSIO_H_

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

/*!
 * Reads a detector list in a .txt file
 * Returns a vector of string containing the name of the considered channels
 */
void read_strings(string fname, std::vector<string>& bolos);

/*!
 * Reads the detectors offsets in a .txt file
 * Returns an array containing the considered channel offsets + the source offsets
 */
void read_bolo_offsets(string field, string file_BoloOffsets, float *scoffsets, double *offsets);

void read_bolo_offsets_from_fits(string filename, string field, double *offsets);

//void read_data_from_fits(string filename, double *data, double *data2, double *data3, short *data4, bool flag, short *data5, long &ns, string field);
//void read_position_from_fits(string filename, double *RA, double *DEC, double *PHI, short *FLAG, bool flag, short *mask, long &ns, string field);
void read_position_from_fits(string filename, double *RA, double *DEC, double *PHI, long &ns);


void read_flpoint_from_fits(string filename, short *FLAG);

void read_flag_from_fits(string filename, short *mask, string field);

void read_signal_from_fits(string filename, double *signal, string field);

#endif /* POSITIONSIO_H_ */
