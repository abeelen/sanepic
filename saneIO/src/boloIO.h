#ifndef BOLOIO_H_
#define BOLOIO_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
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

#endif /* BOLOIO_H_ */
