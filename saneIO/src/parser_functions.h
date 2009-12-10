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

int read_channel_list(dictionary	*ini, std::vector<std::string> &bolonames, int rank);

int read_fits_file_list(dictionary	*ini, struct directories &dir, struct samples &samples_str, int rank);

int read_pixel_size(dictionary	*ini, struct input_commons &com, int rank);

int read_apodize_samples(dictionary	*ini, struct input_commons &com, int rank);

int read_nofillgap(dictionary	*ini, struct input_commons &com, int rank);

int read_box_coord(dictionary	*ini, std::vector<struct box> &boxFile, int rank);

//int read_RA_DEC_min_max(dictionary	*ini, struct user_options_sanepos &u_opt, int &tmpcount);

//int read_RA_DEC_radius_source(dictionary	*ini, struct user_options_sanepos &u_opt, int tmpcount);

int read_map_flagged_data(dictionary	*ini, struct input_commons &com, int rank);

int read_sampling_frequency(dictionary	*ini, struct user_options &u_opt, int rank);

int read_filter_frequency(dictionary	*ini, struct user_options &u_opt, int rank);

int read_noise_cut_freq(dictionary	*ini, std::vector<double> &fcut, int rank);

//int read_noise_file_list(dictionary	*ini, std::vector<string> &extentnoiseSP);

int read_baseline(dictionary	*ini, struct user_options &u_opt, int rank);

int read_correlation(dictionary	*ini, struct user_options &u_opt, int rank);

int read_remove_poly(dictionary	*ini, struct user_options &u_opt, int rank);

int read_projgaps(dictionary	*ini, struct user_options &u_opt, int rank);

int read_iter(dictionary	*ini, int &iterw, int rank);

int read_ell_file(dictionary	*ini, std::string &ellFile, int rank);

int read_map_file(dictionary	*ini, std::string &signame, int rank);

int read_cov_matrix_file(dictionary	*ini, std::string &fname, int rank);

int read_mixmatfile(dictionary	*ini, std::string &MixMatfile, int rank);

int read_directories(dictionary	*ini, struct directories &dir, int rank);

int read_commons(dictionary	*ini, struct input_commons &commons, int rank);

int read_user_options(dictionary *ini,struct user_options &u_opt, int rank);

int read_parser_string(dictionary	*ini, std::string line, int rank, std::string str);

void print_commons(struct input_commons commons);

void print_directories(struct directories dir);

void print_parser(struct user_options u_opt);

#endif /* PARSER_FUNCTIONS_H_ */
