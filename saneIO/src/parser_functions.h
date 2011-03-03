#ifndef PARSER_FUNCTIONS_H_
#define PARSER_FUNCTIONS_H_


#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "utilities.h"
#include "struct_definition.h"

#ifdef USE_MPI
#include "mpi.h"
#endif


extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

std::string checkTrailingDir(std::string str);

void default_param_common(struct param_common &dir);
void default_param_sanePos(struct param_sanePos &pos_param);
void default_param_sanePre(struct param_sanePre &proc_param);
void default_param_sanePS(struct param_sanePS &ps_param);
void default_param_saneInv(struct param_saneInv &inv_param);
void default_param_sanePic(struct param_sanePic &pic_param);

void read_common(std::string &output, dictionary *ini, struct param_common &dir);
void read_param_sanePos(std::string &output, dictionary *ini, struct param_sanePos &pos_param);
void read_param_sanePre(std::string &output, dictionary *ini, struct param_sanePre &proc_param);
void read_param_saneInv(std::string &output, dictionary *ini, struct param_saneInv &saneInv_struct);
void read_param_sanePS(std::string &output,  dictionary *ini, struct param_saneInv &sanePS_struct);
void read_param_sanePic(std::string &output, dictionary *ini, struct param_saneInv &sanePic_struct);

int check_path(std::string &output, std::string strPath, std::string path_type);
int compute_dirfile_format_file(std::string tmp_dir, struct samples samples_struct, int format);
int cleanup_dirfile_sanePos(std::string tmp_dir, struct samples samples_struct);
int cleanup_dirfile_saneInv(std::string tmp_dir, struct samples samples_struct, long n_iter, std::string noise_suffix);
int cleanup_dirfile_fdata(std::string tmp_dir, struct samples samples_struct);

int check_common(std::string &output, struct param_common dir);
int check_param_positions(std::string &output, struct param_sanePos pos_param);
int check_param_process(std::string &output, struct param_sanePre proc_param);
int check_param_sanePS(std::string &output, struct param_sanePS structPS);

void fill_sanePS_struct(struct param_sanePS &structPS, struct samples &samples_struct, struct param_common &dir);
int fill_samples_struct(std::string &output, struct samples &samples_struct, struct param_common &dir, struct param_saneInv &inv_param, std::string fcut_file);
int get_noise_bin_sizes(std::string tmp_dir, struct samples &samples_struct);

#ifdef USE_MPI

void fill_var_sizes_struct(struct param_common dir, struct param_sanePos pos_param, struct param_sanePre proc_param,
		struct param_saneInv inv_param, struct param_sanePS ps_param, struct samples samples_struct, struct ini_var_strings &ini_v);

void Build_derived_type_ini_var (struct ini_var_strings *ini_v,
		MPI_Datatype* message_type_ptr);

int commit_struct_from_root(struct param_common &dir, struct param_sanePos &pos_param, struct param_sanePre &proc_param,
		struct param_saneInv &inv_param, struct param_sanePic &pic_param, struct param_sanePS &ps_param, struct samples &samples_struct, struct ini_var_strings ini_v, int rank);

int commit_param_common(struct param_common &dir, struct ini_var_strings ini_v, int rank);
int commit_param_sanePos(struct param_sanePos &pos_param, struct ini_var_strings ini_v, int rank);
int commit_param_sanePre(struct param_sanePre &proc_param, struct ini_var_strings ini_v, int rank);
int commit_param_saneInv(struct param_saneInv &inv_param, struct ini_var_strings ini_v, int rank);
int commit_param_sanePic(struct param_sanePic &pic_param, struct ini_var_strings ini_v, int rank);
int commit_param_sanePS(struct param_sanePS &ps_param, struct ini_var_strings ini_v, int rank);
int commit_samples_struct(struct samples &samples_struct, struct ini_var_strings ini_v, int rank);


#endif

// this function calls read_* and check_* routines
int parser_function(char * ini_name, std::string &output, struct param_common &dir,
		struct samples &samples_struct,
		struct param_sanePos &pos_param, struct param_sanePre &proc_param,
		struct param_sanePS &structPS, struct param_saneInv &saneInv_struct, struct param_sanePic &sanePic_struct,
		int size, int rank);

void print_common(struct param_common dir);
void print_param_positions(struct param_sanePos pos_param);
void print_param_process(struct param_sanePre proc_param);
void print_param_sanePic(struct param_sanePic sanepic_struct);
void print_param_sanePS(struct param_sanePS structPS);
void print_param_saneInv(struct param_saneInv saneInv_struct);

// calls print_* routines
void parser_printOut(char * prog_name, struct param_common dir, struct samples samples_struct,
		struct param_sanePos pos_param, struct param_sanePre proc_param,
		struct param_sanePS structPS, struct param_sanePic sanePic_struct, struct param_saneInv saneInv_struct);

#endif /* PARSER_FUNCTIONS_H_ */
