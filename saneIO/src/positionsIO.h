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
void read_bolofile(string fname, std::vector<string>& bolos);

/*!
 * Reads the detectors offsets in a .txt file
 * Returns an array containing the considered channel offsets + the source offsets
 */
void read_bolo_offsets(string field, string file_BoloOffsets, float *scoffsets, double *offsets);

void read_bolo_offsets_from_fits(string filename, string field, float *scoffsets, double *offsets);

#endif /* POSITIONSIO_H_ */
