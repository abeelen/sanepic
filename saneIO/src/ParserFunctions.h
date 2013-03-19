#ifndef PARSER_FUNCTIONS_H_
#define PARSER_FUNCTIONS_H_


#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "Utilities.h"
#include "StructDefinition.h"
#include "InputFileIO.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#ifdef USE_MPI
#include "mpi.h"
#endif

//! Given a single value (<= 0 if not used) and a filename, fill a vector<double> of length ntotscan
/*!
  \param value the common positive value
  \param filename a filename containing ntotscan double value to be read
  \param ntotscan the desired output vector size (filename shall contains that exact number of lines...
  \return outputVector the content of filename or value if value > 0
 */

uint32_t fillvect_double(std::string & output, double value, std::string filename, std::string dir, long ntotscan, std::vector<double> &outputVector);

//! Given a common string ("" if not used) and a vector of string and a suffit, fill a vector<string> of length ntotscan
/*!
  \param commonFile if used
  \param Fitsfilename a vector of fits filename
  \param suffix the suffix if commonFile is not used
  \return outputVector contains commonFile if used or FitsBasename(FitsFilename[ii]) + suffix
 */

void fillvect_strings(std::string commonFile, std::vector<std::string> FitsFilename, std::string suffix, std::string dir, std::vector<std::string> & outputVector);

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

void default_param_saneCheck(struct param_saneCheck &Check_param);
void default_param_saneFix(struct param_saneFix &Fix_param);


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

void read_param_saneCheck(std::string &output, dictionary *ini , struct param_saneCheck &Check_param);
void read_param_saneFix(std::string &output, dictionary *ini , struct param_saneFix &Fix_param);

//! Check a path validity or accessibility
/*!
 * Adapt the error message to "path_type" : input, data, output, noise or temporary directory
 * Any warning or error is stored in "output" string and is printed after function exit
 \param strPath A string containing a path name
 \param output The parser error string
 \param path_type A string that determines which is being checked
 \return An integer specifying if there were an error (>0) or not (=0)
 */
uint32_t check_path(std::string &output, std::string strPath, bool create, int rank);

uint32_t init_tmpdir(std::string &output, struct samples &samples_struct, std::string dir, int rank);

//! Check for file existence
/*!
 * \param strPath the complete filepath
 * \return int 0 if file exist
 */
uint32_t check_file(std::string strPath);

//! Check the struct param_common is correct
/*!
 * Any warning or error is stored in "output" string and is printed after function exit
 \param dir The param_common struct
 \param output The parser error string
 \return A flag corresponding to an error code, or 0
 */
uint32_t check_common(std::string &output, struct param_common &dir, int rank);

//! Check the struct param_sanePos is correct
/*!
 * Any warning or error is stored in "output" string and is printed after function exit
 \param pos_param The param_sanePos struct
 \param output The parser error string
 \return A flag corresponding to an error code, or 0
 */
uint32_t check_param_sanePos(std::string &output, struct param_sanePos &pos_param);

//! Check the struct param_saneProc is correct
/*!
 * Any warning or error is stored in "output" string and is printed after function exit
 \param proc_param The param_saneProc struct
 \param output The parser error string
 \return A flag corresponding to an error code, or 0
 */
uint32_t check_param_saneProc(std::string &output, struct param_common dir, struct param_saneProc &proc_param);

//! Check the struct param_sanePS is correct
/*!
 * Any warning or error is stored in "output" string and is printed after function exit
 \param structPS The param_sanePS struct
 \param output The parser error string
 \return A flag corresponding to an error code, or 0
 */
uint32_t check_param_sanePS(std::string &output, struct param_common dir, struct param_sanePS &structPS);

uint32_t check_param_saneInv(std::string &output, struct param_common dir, struct param_saneInv &Inv_param);

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
uint32_t fill_samples_param(std::string &output, struct samples &samples_struct, struct param_common &dir, struct param_saneInv &inv_param, struct param_saneProc &Proc_param, int rank, int size);

//! Open each channel list file and stores those list in a vector (one vector per scan so it's a vector<vector>)
/*!
 \param samples_struct The samples structure that will be filled with nbins and ndet values
 \param bolo_vect The vector containing the channel list of each scan
 \param rank The processor rank given by MPI_Comm_rank, in case paraframe or parabolo is defined
 \return An integer specifying if there were an error (>0) or not (=0)
 */
uint32_t fill_channel_list(std::string &output, struct samples &samples_struct, int rank, int size);

#ifdef USE_MPI

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
uint32_t parser_function(char * ini_name, std::string &output, struct param_common &dir,
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

void print_param_saneCheck(struct param_saneCheck Check_param);
void print_param_saneFix(struct param_saneFix Fix_param);


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
void export_param_saneFix(struct param_saneFix saneFix_struct, std::vector<std::string> &key, std::vector<std::string> &value, std::vector<std::string> &comment);

std::string rebuild_ini(struct param_common dir, struct param_saneProc proc_param, struct param_sanePos pos_param,
		struct samples samples_struct, struct param_sanePS PS_param,
		struct param_sanePic Pic_param, struct param_saneInv Inv_param);

std::string rebuild_ini_saneCheck(struct param_saneCheck Check_param);
std::string rebuild_ini_saneFix(struct param_saneFix Fix_param);

#endif /* PARSER_FUNCTIONS_H_ */
