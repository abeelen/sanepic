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
#include "inputFileIO.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#ifdef USE_MPI
#include "mpi.h"
#endif

//! Parser functions binary return flag
/*! ini file was not found : incorrect name, or name was not given by user */
#define INI_NOT_FOUND 0x0001

//! Parser functions binary return flag
/*! data or/and input path(s) are incorrects or were not correctly filled in ini file */
#define DATA_INPUT_PATHS_PROBLEM 0x0002

//! Parser functions binary return flag
/*! output path is incorrect or was not correctly filled in ini file */
#define OUPUT_PATH_PROBLEM 0x0004

//! Parser functions binary return flag
/*! temporary files path is incorrect or was not correctly filled in ini file */
#define TMP_PATH_PROBLEM 0x0008

//! Parser functions binary return flag
/*! channels files were not found or not correctly filled in ini file */
#define BOLOFILE_NOT_FOUND 0x0010

//! Parser functions binary return flag
/*! pixel size was not correctly filled in ini file or the value is wrong */
#define PIXDEG_WRONG_VALUE 0x0020

//! Parser functions binary return flag
/*! file format was not correctly filled in ini file or the value is absent */
#define FILEFORMAT_NOT_FOUND 0x0040

//! Parser functions binary return flag
/*! napod was not correctly filled in ini file or the value is < 0 */
#define NAPOD_WRONG_VALUE 0x0080

//! Parser functions binary return flag
/*! fsamp was not correctly filled in ini file or the value is < 0 */
#define FSAMP_WRONG_VALUE 0x0100

//! Parser functions binary return flag
/*! f_lp was not correctly filled in ini file or the value is < 0 */
#define F_LP_WRONG_VALUE 0x0200

//! Parser functions binary return flag
/*! ncomp was not correctly filled in ini file or the value is < 0 */
#define NCOMP_WRONG_VALUE 0x0400

//! Parser functions binary return flag
/*! ell file(s) were not found or the name was not correctly given in ini file */
#define ELL_FILE_NOT_FOUND 0x0800

//! Parser functions binary return flag
/*! mixing matrix file(s) were not found or not correctly filled in ini file */
#define MIX_FILE_NOT_FOUND 0x1000

//! Parser functions binary return flag
/*! Noise path or covariance matrices was/were not found or not correctly filled in ini file */
#define SANEINV_INPUT_ERROR 0x2000

//! Parser functions binary return flag
/*! fits filelist file was not found or not correctly filled in ini file */
#define FITS_FILELIST_NOT_FOUND 0x4000

//! Parser functions binary return flag
/*! fcut file was not found or not correctly filled */
#define FCUT_FILE_PROBLEM 0x8000

//! Given a pathname, add a "/", if needed, at the end of the string
/*!
  \param str A path stored in a string
  \return The corrected path (as a string)
 */
std::string checkTrailingDir(std::string str);

//! Fill each struct with default values, in case an optional value was not given in ini file
/*!
 * This routine calls each default_* routine (except sanePS which is not included in the package yet...)
 \param dir The param_common structure
 \param samples_struct The samples structure
 \param pos_param The param_sanePS structure
 \param proc_param The param_saneProc structure
 \param pic_param The param_sanePic structure
 \param inv_param The param_saneInv structure
  \return The default param_common structure
 */
void default_param(struct param_common &dir, struct samples &samples_struct, struct param_sanePos &pos_param, struct param_saneProc &proc_param,
		struct param_saneInv &inv_param, struct param_sanePic &pic_param);


//! Fill the struct param_common with default values, in case an optional value was not given in ini file
/*!
  \param dir An empty param_common structure
  \return The default param_common structure
 */
void default_param_common(struct param_common &dir);


//! Fill the struct param_sanePos with default values, in case an optional value was not given in ini file
/*!
  \param pos_param An empty param_sanePos structure
  \return The default param_sanePos structure
 */
void default_param_sanePos(struct param_sanePos &pos_param);

//! Fill the struct param_saneProc with default values, in case an optional value was not given in ini file
/*!
  \param proc_param An empty param_saneProc structure
  \return The default param_saneProc structure
 */
void default_param_saneProc(struct param_saneProc &proc_param);

//! Fill the struct param_sanePS with default values, in case an optional value was not given in ini file
/*!
  \param ps_param An empty param_sanePS structure
  \return The default param_sanePS structure
 */
void default_param_sanePS(struct param_sanePS &ps_param);

//! Fill the struct param_saneInv with default values, in case an optional value was not given in ini file
/*!
 \param inv_param An empty param_saneInv structure
 \return The default param_saneInv structure
 */
void default_param_saneInv(struct param_saneInv &inv_param);

//! Fill the struct param_sanePic with default values, in case an optional value was not given in ini file
/*!
 \param pic_param An empty param_sanePic structure
 \return The default param_sanePic structure
 */
void default_param_sanePic(struct param_sanePic &pic_param);

//! Fill the struct param_common with ini file values
/*!
 * Any warning or error is stored in "output" string and is printed after function exit
 \param dir The default param_common structure
 \param output The parser error string
 \param ini the ini file opened as a dictionnary with iniparser lib
 \return The filled param_common structure
 */
void read_common(std::string &output, dictionary *ini, struct param_common &dir);


//! Fill the struct param_sanePos with ini file values
/*!
 * Any warning or error is stored in "output" string and is printed after function exit
 \param pos_param The default param_sanePos structure
 \param output The parser error string
 \param ini the ini file opened as a dictionnary with iniparser lib
 \return The filled param_sanePos structure
 */
void read_param_sanePos(std::string &output, dictionary *ini, struct param_sanePos &pos_param);

//! Fill the struct param_saneProc with ini file values
/*!
 * Any warning or error is stored in "output" string and is printed after function exit
 \param proc_param The default param_saneProc structure
 \param output The parser error string
 \param ini the ini file opened as a dictionnary with iniparser lib
 \return The filled param_saneProc structure
 */
void read_param_saneProc(std::string &output, dictionary *ini, struct param_saneProc &proc_param);

//! Fill the struct param_saneInv with ini file values
/*!
 * Any warning or error is stored in "output" string and is printed after function exit
 \param saneInv_struct The default param_saneInv structure
 \param output The parser error string
 \param ini the ini file opened as a dictionnary with iniparser lib
 \return The filled param_saneInv structure
 */
void read_param_saneInv(std::string &output, dictionary *ini, struct param_saneInv &saneInv_struct);

//! Fill the struct param_sanePS with ini file values
/*!
 * Any warning or error is stored in "output" string and is printed after function exit
 \param sanePS_struct The default param_sanePS structure
 \param output The parser error string
 \param ini the ini file opened as a dictionnary with iniparser lib
 \return The filled param_sanePS structure
 */
void read_param_sanePS(std::string &output,  dictionary *ini, struct param_sanePS &sanePS_struct);

//! Fill the struct param_sanePic with ini file values
/*!
 * Any warning or error is stored in "output" string and is printed after function exit
 \param sanePic_struct The default param_sanePic structure
 \param output The parser error string
 \param ini the ini file opened as a dictionnary with iniparser lib
 \return The filled param_sanePic structure
 */
void read_param_sanePic(std::string &output, dictionary *ini, struct param_sanePic &sanePic_struct);

//! Check a path validity or accessibility
/*!
 * Adapt the error message to "path_type" : input, data, output, noise or temporary directory
 * Any warning or error is stored in "output" string and is printed after function exit
 \param strPath A string containing a path name
 \param output The parser error string
 \param path_type A string that determines which is being checked
 \return An integer specifying if there were an error (>0) or not (=0)
 */
int check_path(std::string &output, std::string strPath, std::string path_type, bool create);

//! Creates dirfile architecture and format files considering fits file format (SANEPIC or HIPE)
/*!
 * Each scan has his own branch, with : data / flag / Indexes / fData / Noise_data
 \param format An integer : HIPE (1) or SANEPIC (2)
 \param samples_struct A samples structure
 \param tmp_dir A string containing the temporary files pathname
 \return An integer specifying if there were an error (>0) or not (=0)
 */
int compute_dirfile_format_file(std::string tmp_dir, struct samples samples_struct, int format);

//! Clean up Indexes dirfiles and format files or Creates it if needed
/*!
 * Each scan has his own Indexes dirfile
 \param bolo_vect A vector containing the channel list (as a vector of string), for whole scan
 \param samples_struct A samples structure
 \param tmp_dir A string containing the temporary files pathname
 \return An integer specifying if there were an error (>0) or not (=0)
 */
int cleanup_dirfile_sanePos(std::string tmp_dir, struct samples samples_struct, std::vector<std::vector<std::string> > bolo_vect);

//! Clean up Noise_data and Noise_data/ell dirfiles and format files or Creates it if needed
/*!
 * Each scan has his own Noise_data dirfile
 \param bolo_vect A vector containing the channel list (as a vector of string), for whole scan
 \param noise_suffix A suffix to add to each fits file Basename to form the noise filename
 \param samples_struct A samples structure
 \param nframe The number of scans
 \param tmp_dir A string containing the temporary files pathname
 \return An integer specifying if there were an error (>0) or not (=0)
 */
int cleanup_dirfile_saneInv(std::string tmp_dir, struct samples samples_struct, long nframe, std::string noise_suffix, std::vector<std::vector<std::string> > bolo_vect);

//! Clean up fData dirfiles and format files or Creates it if needed
/*!
 * Each scan has his own fData dirfile
 \param bolo_vect A vector containing the channel list (as a vector of string), for whole scan
 \param samples_struct A samples structure
 \param tmp_dir A string containing the temporary files pathname
 \return An integer specifying if there were an error (>0) or not (=0)
 */
int cleanup_dirfile_fdata(std::string tmp_dir, struct samples samples_struct, std::vector<std::vector<std::string> > bolo_vect);

//! Check the struct param_common is correct
/*!
 * Any warning or error is stored in "output" string and is printed after function exit
 \param dir The param_common struct
 \param output The parser error string
 \return A flag corresponding to an error code, or 0
 */
uint16_t check_common(std::string &output, struct param_common dir);

//! Check the struct param_sanePos is correct
/*!
 * Any warning or error is stored in "output" string and is printed after function exit
 \param pos_param The param_sanePos struct
 \param output The parser error string
 \return A flag corresponding to an error code, or 0
 */
uint16_t check_param_sanePos(std::string &output, struct param_sanePos pos_param);

//! Check the struct param_saneProc is correct
/*!
 * Any warning or error is stored in "output" string and is printed after function exit
 \param proc_param The param_saneProc struct
 \param output The parser error string
 \return A flag corresponding to an error code, or 0
 */
uint16_t check_param_saneProc(std::string &output, struct param_saneProc proc_param);

//! Check the struct param_sanePS is correct
/*!
 * Any warning or error is stored in "output" string and is printed after function exit
 \param structPS The param_sanePS struct
 \param output The parser error string
 \return A flag corresponding to an error code, or 0
 */
uint16_t check_param_sanePS(std::string &output, struct param_sanePS structPS);

//! Fill param_sanePS structure with ell and mixing matrices names
/*!
 \param structPS The param_sanePS struct partially filled to be completed
 \param samples_struct The samples structure that will be used to fill structPS
 \param dir The param_common struct
 \return The filled param_sanePS struct
 */
void fill_sanePS_struct(struct param_sanePS &structPS, struct samples &samples_struct, struct param_common &dir);

//! Fill samples structure with ini, fcut and fitsfilelist files informations
/*!
 \param output The parser error string
 \param dir The param_common structure
 \param samples_struct The samples structure that will be filled
 \param inv_param The param_saneInv struct
 \param fcut_file A string containing the fcut filename
 \return A flag corresponding to an error code, or 0
 */
uint16_t fill_samples_struct(std::string &output, struct samples &samples_struct, struct param_common &dir, struct param_saneInv &inv_param, std::string fcut_file);

//! Open a dirfile to get nbins (Noise_data/ell) and ndet (Noise_data) values
/*!
 \param tmp_dir A string containing the temporary files path
 \param samples_struct The samples structure that will be filled with nbins and ndet values
 \param rank The processor rank given by MPI_Comm_rank, in case paraframe or parabolo is defined
 \return An integer specifying if there were an error (>0) or not (=0)
 */
int get_noise_bin_sizes(std::string tmp_dir, struct samples &samples_struct, int rank);

//! Open each channel list file and stores those list in a vector (one vector per scan so it's a vector<vector>)
/*!
 \param samples_struct The samples structure that will be filled with nbins and ndet values
 \param bolo_vect The vector containing the channel list of each scan
 \param rank The processor rank given by MPI_Comm_rank, in case paraframe or parabolo is defined
 \return An integer specifying if there were an error (>0) or not (=0)
 */
int channel_list_to_vect_list(struct samples samples_struct, std::vector<std::vector<std::string> > &bolo_vect, int rank);

#ifdef USE_MPI

//! Computes the size of each string contained in "str_vect"
/*!
 * Stores the maximum size in size_max and returns the whole strings sizes after summation
 \param str_vect A vector containing a list of channels (their names)
 \param size_max A long int that corresponds to the longest channel's name
 \return size_buff : A long int that corresponds to the summation of whole string sizes
 */
long compute_bololist_size(std::vector<std::string> str_vect, long &size_max);


int commit_dictionary(int rank, dictionary	*dict);

#endif

//! Fill structures with ini, fcut and fitsfilelist files informations
/*!
 * This function calls default_*, read_* and check_* routines
 \param ini_name the ini file name given by user in the command line (usually argv[1])
 \param output The parser error string
 \param dir The param_common structure
 \param samples_struct The samples structure
 \param pos_param The param_sanePos structure
 \param proc_param The param_saneProc structure
 \param structPS The param_sanePS structure
 \param sanePic_struct The param_sanePic structure
 \param saneInv_struct The param_saneInv structure
 \param size Total number of processors, in case paraframe or parabolo is defined
 \param rank The processor rank given by MPI_Comm_rank, in case paraframe or parabolo is defined
 \return A flag corresponding to an error code, or 0
 */
uint16_t parser_function(char * ini_name, std::string &output, struct param_common &dir,
		struct samples &samples_struct,
		struct param_sanePos &pos_param, struct param_saneProc &proc_param,
		struct param_sanePS &structPS, struct param_saneInv &saneInv_struct, struct param_sanePic &sanePic_struct,
		int size, int rank);

//! Print param_common structure informations to screen
/*!
 \param dir The param_common structure
 */
void print_common(struct param_common dir);

//! Print param_sanePos structure informations to screen
/*!
 \param pos_param The param_sanePos structure
 */
void print_param_sanePos(struct param_sanePos pos_param);

//! Print param_saneProc structure informations to screen
/*!
 \param proc_param The param_saneProc structure
 */
void print_param_saneProc(struct param_saneProc proc_param);

//! Print param_sanePic structure informations to screen
/*!
 \param sanepic_struct The param_sanePic structure
 */
void print_param_sanePic(struct param_sanePic sanepic_struct);

//! Print param_sanePS structure informations to screen
/*!
 \param structPS The param_sanePS structure
 */
void print_param_sanePS(struct param_sanePS structPS);

//! Print param_saneInv structure informations to screen
/*!
 \param saneInv_struct The param_saneInv structure
 */
void print_param_saneInv(struct param_saneInv saneInv_struct);

//! Print structures informations to screen depending on which program is calling it
/*!
 * Call the functions print_*
 \param prog_name The program that is currently calling the function : argv[0]
 \param dir The param_common structure
 \param samples_struct The samples structure
 \param pos_param The param_sanePos structure
 \param proc_param The param_saneProc structure
 \param structPS The param_sanePS structure
 \param sanePic_struct The param_sanePic structure
 \param saneInv_struct The param_saneInv structure
 */
void parser_printOut(char * prog_name, struct param_common dir, struct samples samples_struct,
		struct param_sanePos pos_param, struct param_saneProc proc_param,
		struct param_sanePS structPS, struct param_sanePic sanePic_struct, struct param_saneInv saneInv_struct);



void export_param_sanePos(struct param_sanePos pos_param, std::vector<std::string> &key, std::vector<std::string> &value, std::vector<std::string> &comment);

void export_param_common(struct param_common dir, std::vector<std::string> &key, std::vector<std::string> &value, std::vector<std::string> &comment);

void export_param_saneProc(struct param_saneProc proc_param, std::vector<std::string> &key, std::vector<std::string> &value, std::vector<std::string> &comment);

void export_param_saneInv(struct param_saneInv inv_param, std::vector<std::string> &key, std::vector<std::string> &value, std::vector<std::string> &comment);

void export_param_sanePS(struct param_sanePS ps_param, std::vector<std::string> &key, std::vector<std::string> &value, std::vector<std::string> &comment);

void export_param_sanePic(struct param_sanePic sanePic_struct, std::vector<std::string> &key, std::vector<std::string> &value, std::vector<std::string> &comment);

void export_param_saneCheck(struct param_saneCheck saneCheck_struct, std::vector<std::string> &key, std::vector<std::string> &value, std::vector<std::string> &comment);

std::string rebuild_ini(struct param_common dir, struct param_saneProc proc_param, struct param_sanePos pos_param,
		struct samples samples_struct, struct param_sanePS PS_param,
		struct param_sanePic Pic_param, struct param_saneInv Inv_param);

#endif /* PARSER_FUNCTIONS_H_ */
