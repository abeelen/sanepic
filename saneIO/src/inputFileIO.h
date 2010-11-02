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

std::string FitsBasename(std::string path);

int read_bolo_for_all_scans(std::vector<detectors> &detector_tab, struct param_common dir, struct samples samples_struct, int rank, int size);
int read_channel_list(std::string &output, std::string fname, std::vector<std::string> &bolonames);

#endif /* INPUTFILEIO_H_ */
