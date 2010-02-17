/*
 * parser_functions.h
 *
 *  Created on: 20 oct. 2009
 *      Author: matthieu
 */

#ifndef PARSER_FUNCTIONS_H_
#define PARSER_FUNCTIONS_H_


#include <string>
#include <vector>
#include "struct_definition.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}



int read_dirfile(dictionary	*ini, struct directories &dir, int rank);
int read_tmpdir(dictionary	*ini, struct directories &dir, int rank);
int read_outdir(dictionary	*ini, struct directories &dir, int rank);
int read_outdir(dictionary	*ini, struct directories &dir, int rank);

int read_channel_list(dictionary	*ini, std::vector<std::string> &bolonames, int rank);
int read_fits_file_list(dictionary	*ini, struct directories &dir, struct samples &samples_str, int rank);
int read_fits_list(std::string fname, std::vector<std::string> &fitsfiles, std::vector<std::string> &noisefiles, std::vector<int> &frameorder, bool &framegiven);

int read_apodize_samples(dictionary	*ini, struct param_process &com, int rank);
int read_nofillgap(dictionary	*ini, struct param_process &com, int rank);

int read_sampling_frequency(dictionary	*ini, struct param_process &proc_param, int rank);

int read_filter_frequency(dictionary	*ini, struct param_process &proc_param, int rank);

int read_noise_cut_freq(dictionary	*ini, std::vector<double> &fcut, int rank);

//int read_noise_file_list(dictionary	*ini, std::vector<string> &extentnoiseSP);
int read_baseline(dictionary	*ini, struct param_process &proc_param, int rank);
int read_correlation(dictionary	*ini, struct param_process &proc_param, int rank);
int read_remove_poly(dictionary	*ini, struct param_process &proc_param, int rank);
int read_iter(dictionary	*ini, int &iterw, int rank);
int read_ell_file(dictionary	*ini, std::string &ellFile, int rank);
int read_map_file(dictionary	*ini, std::string &signame, int rank);
int read_cov_matrix_file(dictionary	*ini, std::string &fname, int rank);
int read_mixmatfile(dictionary	*ini, std::string &MixMatfile, int rank);
int read_ncomp(dictionary	*ini, long &ncomp, int rank);
int read_fcut(dictionary	*ini, double &fcut, int rank);

int read_directories(dictionary	*ini, struct directories &dir, int rank);

int read_parser_string(dictionary	*ini, std::string line, int rank, std::string & str);
int read_param_process(dictionary *ini,struct param_process &proc_param, int rank);
int read_param_positions(dictionary *ini, struct param_positions &pos_param, int rank);

void print_param_positions(struct param_positions pos_param);
void print_param_process(struct param_process proc_param);

void print_directories(struct directories dir);

#endif /* PARSER_FUNCTIONS_H_ */
