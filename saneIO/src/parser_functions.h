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
int read_parser_string(dictionary *ini,  std::string line, std::string & str);
int read_fits_file_list(std::string &output, dictionary	*ini, struct param_common &dir, struct samples &samples_str);

int read_noise_cut_freq(std::string &output, dictionary	*ini, struct param_common dir, struct param_sanePre &proc_param, std::vector<double> &fcut);
int read_ell_suffix(std::string &output, dictionary	*ini, std::string &ell_suffix);
int read_ell_global_file(std::string &output, dictionary	*ini, std::string &ell_global_file);
int read_map_file(dictionary *ini, std::string &signame);
int read_cov_matrix_file(std::string &output, dictionary *ini, std::string &fname);
int read_cov_matrix_suffix(std::string &output, dictionary *ini, std::string &fname);
int read_mixmatfile_suffix(std::string &output, dictionary *ini, std::string &MixMat_suffix);
int read_mixmat_global_file(std::string &output, dictionary	*ini, std::string &MixMat_global);
int read_bolo_suffix(dictionary	*ini, std::string &suffix);
int read_bolo_global_file(dictionary *ini, std::string &bolo_global_filename);
int read_bolo_gain_global_file(std::string &output, dictionary *ini, std::string dir, std::string &bolo_global_filename);

int read_common(std::string &output, dictionary	*ini, struct param_common &dir);
int read_param_positions(std::string &output, dictionary *ini, struct param_common dir, struct param_sanePos &pos_param);
int read_param_process(std::string &output, dictionary *ini, struct param_common dir, struct param_sanePre &proc_param, std::vector<double> &fcut);
int read_param_saneInv(std::string &output, dictionary *ini, struct param_saneInv &saneInv_struct);
int read_param_sanePS(std::string &output, dictionary *ini, struct param_saneInv &sanePS_struct);
int read_param_sanePic(std::string &output, dictionary *ini, struct param_saneInv &sanePic_struct);

int check_path(std::string &output, std::string strPath, std::string path_type);
int check_dirfile_paths(std::string &output, std::string strPath);

int check_common(std::string &output, struct param_common dir);
int check_param_positions(std::string &output, struct param_sanePos pos_param);
int check_param_process(std::string &output, struct param_sanePre proc_param);
int check_param_sanePS(std::string &output, struct param_sanePS structPS);

// this function calls read_* and check_* routines
int parser_function(char * ini_name, std::string &output, struct param_common &dir,
		struct samples &samples_struct,
		struct param_sanePos &pos_param, struct param_sanePre &proc_param, std::vector<double> &fcut,
		struct param_sanePS &structPS, struct param_saneInv &saneInv_struct, struct param_sanePic &sanePic_struct,
		int size);

void print_common(struct param_common dir);
void print_param_positions(struct param_sanePos pos_param);
void print_param_process(struct param_sanePre proc_param);
void print_param_sanePic(struct param_sanePic sanepic_struct);
void print_param_sanePS(struct param_sanePS structPS);

// calls print_* routines
void parser_printOut(char * prog_name, struct param_common dir, struct samples samples_struct,
		struct param_sanePos pos_param, struct param_sanePre proc_param,
		struct param_sanePS structPS, struct param_sanePic sanePic_struct);

#endif /* PARSER_FUNCTIONS_H_ */
