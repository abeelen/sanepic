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


extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#ifdef USE_MPI
#include "mpi.h"
#endif

#define INI_NOT_FOUND 0x0001
#define DATA_INPUT_PATHS_PROBLEM 0x0002
#define OUPUT_PATH_PROBLEM 0x0004
#define TMP_PATH_PROBLEM 0x0008
#define BOLOFILE_NOT_FOUND 0x0010
#define PIXDEG_WRONG_VALUE 0x0020
#define FILEFORMAT_NOT_FOUND 0x0040
#define NAPOD_WRONG_VALUE 0x0080
#define FSAMP_WRONG_VALUE 0x0100
#define F_LP_WRONG_VALUE 0x0200
#define NCOMP_WRONG_VALUE 0x0400
#define ELL_FILE_NOT_FOUND 0x0800
#define MIX_FILE_NOT_FOUND 0x1000
// noise dir or cov matrix not found
#define SANEINV_INPUT_ERROR 0x2000
#define FITS_FILELIST_NOT_FOUND 0x4000
// file not found or fcut values error
#define FCUT_FILE_PROBLEM 0x8000



//ini_not_found + data_input_path_problem + output_path_problem + tmp_path_problem +  bolofile_not_found +
//pixdeg_wrong_value + fileformat_not_found + napod_wrong_value + fsamp_wrong_value + f_lp_wrong_value +
//ncomp_wrong_value + ell_file_not_found + mix_file_not_found + saneInv_input_error + fits_filelist_not_found+
//fcut_file_problem;

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
int cleanup_dirfile_sanePos(std::string tmp_dir, struct samples samples_struct, std::vector<std::vector<std::string> > bolo_vect);
int cleanup_dirfile_saneInv(std::string tmp_dir, struct samples samples_struct, long n_iter, std::string noise_suffix, std::vector<std::vector<std::string> > bolo_vect);
int cleanup_dirfile_fdata(std::string tmp_dir, struct samples samples_struct, std::vector<std::vector<std::string> > bolo_vect);

uint16_t check_common(std::string &output, struct param_common dir);
uint16_t check_param_positions(std::string &output, struct param_sanePos pos_param);
uint16_t check_param_process(std::string &output, struct param_sanePre proc_param);
uint16_t check_param_sanePS(std::string &output, struct param_sanePS structPS);

void fill_sanePS_struct(struct param_sanePS &structPS, struct samples &samples_struct, struct param_common &dir);
uint16_t fill_samples_struct(std::string &output, struct samples &samples_struct, struct param_common &dir, struct param_saneInv &inv_param, std::string fcut_file);
int get_noise_bin_sizes(std::string tmp_dir, struct samples &samples_struct, int rank);

//int channel_list_to_chain_list(struct samples samples_struct, struct bolo_chaine *ptr, int rank);
int channel_list_to_vect_list(struct samples samples_struct, std::vector<std::vector<std::string> > &bolo_vect, int rank);
long compute_bololist_size(std::vector<std::string> str_vect, long &size_max);

#ifdef USE_MPI

void fill_var_sizes_struct(struct param_common dir, struct param_sanePos pos_param, struct param_sanePre proc_param,
		struct param_saneInv inv_param, struct param_sanePS ps_param, struct samples samples_struct, struct ini_var_strings &ini_v);

void Build_derived_type_ini_var (struct ini_var_strings *ini_v,
		MPI_Datatype* message_type_ptr);

//long compute_bololist_size(std::string **bololist, long iframe, long &size_max, long *size_tab, long *ndet_tab);
//int commit_bololist(std::string **bololist, long iframe, long size_buff, long size_max, long ndet, long *size_tab, int rank);
//int commit_bololist(std::string **bololist, long ntotscan, long *ndet_tab, int rank);
//long compute_bololist_size(std::vector<std::string> str_vect, long &size_max);

//int commit_bololist(std::vector<std::string> &str_vect, long ntotscan, long nbolo, int rank);

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
uint16_t parser_function(char * ini_name, std::string &output, struct param_common &dir,
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
