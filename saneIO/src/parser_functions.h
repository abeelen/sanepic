#ifndef PARSER_FUNCTIONS_H_
#define PARSER_FUNCTIONS_H_


#include <string>
#include <vector>
#include "struct_definition.h"


extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

//template <typename T>
//std::string StringOf(const T& object);

int read_dir(std::string &output, dictionary	*ini, struct param_common &dir, std::string dirtype ,int rank);


int read_fits_file_list(std::string &output, dictionary	*ini, struct param_common &dir, struct samples &samples_str, int rank);
int read_fits_list(std::string &output, std::string fname, std::vector<std::string> &fitsfiles, std::vector<int> &frameorder, bool &framegiven);

int read_apodize_samples(std::string &output, dictionary	*ini, struct param_sanePre &com, int rank);

int read_sampling_frequency(std::string &output, dictionary	*ini, struct param_sanePre &proc_param, int rank);

int read_filter_frequency(std::string &output, dictionary	*ini, struct param_sanePre &proc_param, int rank);

int read_fcut(std::string &output, dictionary	*ini, double &fcut, int rank);
int read_noise_cut_freq(std::string &output, dictionary	*ini, struct param_common dir, struct param_sanePre &proc_param, std::vector<double> &fcut, int rank);
int read_baseline(dictionary	    *ini, struct param_sanePre &proc_param, int rank);
int read_correlation(dictionary	    *ini, struct param_sanePre &proc_param, int rank);
int read_remove_poly(dictionary	    *ini, struct param_sanePre &proc_param, int rank);
int read_iter(dictionary	        *ini, int &iterw, int rank);
int read_ell_suffix(std::string &output, dictionary	*ini, std::string &ell_suffix, int rank);
int read_ell_global_file(std::string &output, dictionary	*ini, std::string &ell_global_file, int rank);
int read_map_file(dictionary	*ini, std::string &signame);
int read_cov_matrix_file(std::string &output, dictionary	*ini, std::string &fname, int rank);
int read_cov_matrix_suffix(std::string &output, dictionary	*ini, std::string &fname, int rank);
int read_mixmatfile_suffix(dictionary	*ini, std::string &MixMat_suffix, int rank);
int read_mixmat_global_file(dictionary	*ini, std::string &MixMat_global, int rank);
int read_ncomp(std::string &output, dictionary	*ini, long &ncomp, int rank);

int read_nofillgap(dictionary	*ini, struct param_sanePre &com, int rank);

int read_common(std::string &output, dictionary	*ini, struct param_common &dir, int rank);

int read_parser_string(dictionary	*ini,  std::string line, std::string & str);
int read_param_process(std::string &output, dictionary *ini,struct param_sanePre &proc_param, int rank);
int read_param_positions(std::string &output, dictionary *ini, struct param_sanePos &pos_param, int rank);

void print_param_positions(struct param_sanePos pos_param);
void print_param_process(struct param_sanePre proc_param);

void print_common(struct param_common dir);

int check_path(std::string &output, std::string strPath, std::string path_type);
int check_dirfile_paths(std::string &output, std::string strPath);

int read_bolo_suffix(dictionary	*ini, std::string &suffix);
int read_bolo_global_file(dictionary *ini, std::string &bolo_global_filename);
int read_bolo_gain_global_file(std::string &output, dictionary *ini, std::string dir, std::string &bolo_global_filename, int rank);

void fill_sanePS_struct(struct param_sanePS &structPS, struct samples samples_struct);
void fill_noisevect(struct samples &samples_struct);

int parser_function(char * ini_name, std::string &output, struct param_common &dir,
		struct samples &samples_struct,
		struct param_sanePos &pos_param, struct param_sanePre &proc_param, std::vector<double> &fcut,
		struct param_sanePS &structPS, struct param_sanePic &sanePic_struct, int rank, int size);

void parser_printOut(struct param_common dir, struct samples samples_struct,
		struct param_sanePos pos_param, struct param_sanePre proc_param,
		struct param_sanePS structPS, struct param_sanePic sanePic_struct, int rank);

#endif /* PARSER_FUNCTIONS_H_ */
