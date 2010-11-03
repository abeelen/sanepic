#ifndef PARSER_FUNCTIONS_H_
#define PARSER_FUNCTIONS_H_


#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "struct_definition.h"


extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

template <typename T>
std::string StringOf(const T& object){
	std::ostringstream os;
	os << object;
	return os.str();
}

int read_dir(std::string &output, dictionary	*ini, struct param_common &dir, std::string dirtype);


int read_fits_file_list(std::string &output, dictionary	*ini, struct param_common &dir, struct samples &samples_str);

int read_apodize_samples(std::string &output, dictionary	*ini, struct param_sanePre &com);

int read_sampling_frequency(std::string &output, dictionary	*ini, struct param_sanePre &proc_param);

int read_filter_frequency(std::string &output, dictionary	*ini, struct param_sanePre &proc_param);

int read_fcut(std::string &output, dictionary	*ini, double &fcut);
int read_noise_cut_freq(std::string &output, dictionary	*ini, struct param_common dir, struct param_sanePre &proc_param, std::vector<double> &fcut);
int read_baseline(std::string &output, dictionary *ini, struct param_sanePre &proc_param);
int read_correlation(std::string &output, dictionary	*ini, struct param_sanePre &proc_param);
int read_remove_poly(std::string &output, dictionary	*ini, struct param_sanePre &proc_param);
int read_iter(std::string &output, dictionary *ini, int &iterw);
int read_ell_suffix(std::string &output, dictionary	*ini, std::string &ell_suffix);
int read_ell_global_file(std::string &output, dictionary	*ini, std::string &ell_global_file);
int read_map_file(dictionary *ini, std::string &signame);
int read_cov_matrix_file(std::string &output, dictionary *ini, std::string &fname);
int read_cov_matrix_suffix(std::string &output, dictionary *ini, std::string &fname);
int read_mixmatfile_suffix(dictionary *ini, std::string &MixMat_suffix);
int read_mixmat_global_file(dictionary	*ini, std::string &MixMat_global);
int read_ncomp(std::string &output, dictionary	*ini, long &ncomp);

int read_nofillgap(std::string &output, dictionary *ini, struct param_sanePre &com);

int read_common(std::string &output, dictionary	*ini, struct param_common &dir);

int read_parser_string(dictionary *ini,  std::string line, std::string & str);
int read_param_process(std::string &output, dictionary *ini, struct param_common dir, struct param_sanePre &proc_param, std::vector<double> &fcut);
int read_param_positions(std::string &output, dictionary *ini, struct param_sanePos &pos_param);
int read_param_saneInv(std::string &output, dictionary *ini, struct param_saneInv &saneInv_struct);
int read_param_sanePic(std::string &output, dictionary *ini, struct param_saneInv &sanePic_struct);
int read_param_sanePS(std::string &output, dictionary *ini, struct param_saneInv &sanePS_struct);

void print_param_positions(struct param_sanePos pos_param);
void print_param_process(struct param_sanePre proc_param);

void print_common(struct param_common dir);

int check_path(std::string &output, std::string strPath, std::string path_type);
int check_dirfile_paths(std::string &output, std::string strPath);

int read_bolo_suffix(dictionary	*ini, std::string &suffix);
int read_bolo_global_file(dictionary *ini, std::string &bolo_global_filename);
int read_bolo_gain_global_file(std::string &output, dictionary *ini, std::string dir, std::string &bolo_global_filename);

void fill_sanePS_struct(struct param_sanePS &structPS, struct samples samples_struct);
void fill_noisevect_fcut(struct param_common dir, struct samples &samples_str, struct param_saneInv &saneInv_struct, std::vector<double> &fcut);

int parser_function(char * ini_name, std::string &output, struct param_common &dir,
		struct samples &samples_struct,
		struct param_sanePos &pos_param, struct param_sanePre &proc_param, std::vector<double> &fcut,
		struct param_sanePS &structPS, struct param_saneInv &saneInv_struct, struct param_sanePic &sanePic_struct,
		int size);

void parser_printOut(struct param_common dir, struct samples samples_struct, std::vector<detectors> detector_tab,
		struct param_sanePos pos_param, struct param_sanePre proc_param,
		struct param_sanePS structPS, struct param_sanePic sanePic_struct);

#endif /* PARSER_FUNCTIONS_H_ */
