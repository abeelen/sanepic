#ifndef PARSER_FUNCTIONS_H_
#define PARSER_FUNCTIONS_H_


#include <string>
#include <vector>
#include "struct_definition.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}


int read_dir(dictionary	*ini, struct common &dir, std::string dirtype ,int rank);

int read_channel_list(std::string fname, std::vector<std::string> &bolonames, int rank);
int read_fits_file_list(dictionary	*ini, struct common &dir, struct samples &samples_str, int rank);
int read_fits_list(std::string fname, std::vector<std::string> &fitsfiles, std::vector<std::string> &noisefiles, std::vector<int> &frameorder, bool &framegiven);

int read_apodize_samples(dictionary	*ini, struct param_process &com, int rank);

int read_sampling_frequency(dictionary	*ini, struct param_process &proc_param, int rank);

int read_filter_frequency(dictionary	*ini, struct param_process &proc_param, int rank);


int read_fcut(dictionary	*ini, double &fcut, int rank);
int read_noise_cut_freq(dictionary	*ini, struct common dir, struct param_process &proc_param, std::vector<double> &fcut, int rank);
int read_baseline(dictionary	    *ini, struct param_process &proc_param, int rank);
int read_correlation(dictionary	    *ini, struct param_process &proc_param, int rank);
int read_remove_poly(dictionary	    *ini, struct param_process &proc_param, int rank);
int read_iter(dictionary	        *ini, int &iterw, int rank);
int read_ell_suffix(dictionary	*ini, std::string &ell_suffix, int rank);
int read_ell_global_file(dictionary	*ini, std::string &ell_global_file, int rank);
int read_map_file(dictionary	*ini, std::string &signame);
int read_cov_matrix_file(dictionary	*ini, std::string &fname, int rank);
int read_mixmatfile_suffix(dictionary	*ini, std::string &MixMat_suffix, int rank);
int read_mixmat_global_file(dictionary	*ini, std::string &MixMat_global, int rank);
int read_ncomp(dictionary	*ini, long &ncomp, int rank);

int read_nofillgap(dictionary	*ini, struct param_process &com, int rank);


int read_common(dictionary	*ini, struct common &dir, int rank);

int read_parser_string(dictionary	*ini,  std::string line, std::string & str);
int read_param_process(dictionary *ini,struct param_process &proc_param, int rank);
int read_param_positions(dictionary *ini, struct param_positions &pos_param, int rank);

void print_param_positions(struct param_positions pos_param);
void print_param_process(struct param_process proc_param);

void print_common(struct common dir);

int check_path(std::string strPath, std::string path_type);
int check_dirfile_paths(std::string strPath);

int read_bolo_suffix(dictionary	*ini, std::string &suffix);
int read_bolo_global_file(dictionary *ini, std::string &bolo_global_filename);
int read_bolo_gain_global_file(dictionary *ini, std::string dir, std::string &bolo_global_filename, int rank);

void fill_sanePS_struct(std::string dir, struct PS &structPS, struct samples samples_struct);

int parser_function(char * ini_name, struct common &dir,
		std::vector<detectors> &detector_tab,struct samples &samples_struct,
		struct param_positions &pos_param, struct param_process &proc_param, std::vector<double> &fcut,
		struct PS &structPS, struct sanePic &sanePic_struct, int rank, int size);

#endif /* PARSER_FUNCTIONS_H_ */
