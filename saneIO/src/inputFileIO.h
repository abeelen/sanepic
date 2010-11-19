#ifndef INPUTFILEIO_H_
#define INPUTFILEIO_H_

#include <string>
#include <vector>
#include "struct_definition.h"

/*!
 * Reads a detector list in a .txt file
 * Returns a vector of string containing the name of the considered channels
 */
int read_strings(std::string fname, std::vector<std::string>& bolos);
int read_double(std::string fname, double *& array, long & size);

std::string remplace_all(std::string str, std::string tobe_replace, std::string with_this);
std::string FitsBasename(std::string path);

void readFrames(std::vector<std::string> &inputFiles, long *&nsamples);

long readFitsLength(std::string filename);
int read_channel_list(std::string &output, std::string fname, std::vector<std::string> &bolonames);
int read_fits_list(std::string &output, std::string fname, struct samples &samples_struct);



#endif /* INPUTFILEIO_H_ */
