/*
 * inputFileIO.h
 */

#ifndef INPUTFILEIO_H_
#define INPUTFILEIO_H_

#include <string>
#include <vector>

using namespace std;

/*!
 * Reads a detector list in a .txt file
 * Returns a vector of string containing the name of the considered channels
 */
void read_strings(string fname, std::vector<string>& bolos);
void read_double(string fname, double *& array, long & size);

#endif /* INPUTFILEIO_H_ */
