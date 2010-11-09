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
//int read_bolo_for_all_scans(std::vector<detectors> &detector_tab, struct param_common dir, struct samples samples_struct, int rank, int size);
int read_bolo_for_all_scans(struct param_common dir, struct samples &samples_struct, int rank, int size);
int read_channel_list(std::string &output, std::string fname, std::vector<std::string> &bolonames);
int read_fits_list(std::string &output, std::string fname, std::vector<std::string> &fitsfiles, std::vector<int> &frameorder, bool &framegiven);

void fill_sanePS_struct(struct param_sanePS &structPS, struct samples samples_struct);
void fill_noisevect_fcut(struct param_common dir, struct samples &samples_str, struct param_saneInv &saneInv_struct, std::vector<double> &fcut);


#endif /* INPUTFILEIO_H_ */
