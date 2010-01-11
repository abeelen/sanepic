/*
 * inputFileIO.h
 */

#ifndef INPUTFILEIO_H_
#define INPUTFILEIO_H_

#include <string>
#include <vector>


/*!
 * Reads a detector list in a .txt file
 * Returns a vector of string containing the name of the considered channels
 */
void read_strings(std::string fname, std::vector<std::string>& bolos);
void read_double(std::string fname, double *& array, long & size);

std::string		Basename(std::string path);
std::string FitsBasename(std::string path);

#endif /* INPUTFILEIO_H_ */
