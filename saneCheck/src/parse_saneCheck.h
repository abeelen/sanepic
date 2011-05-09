#include "struct_definition.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#ifndef PARSE_SANECHECK_H_

//! Fill structures with ini, fcut and fitsfilelist files informations. Generates a ini file model in case some lines were missing in input one
/*!
 * This function calls parser_function(...) (from saneIO), read_saneCheck_ini(...) and print_saneCheck_ini(...) routines
 \param ini_name the ini file name given by user in the command line (usually argv[1])
 \param output The parser error string
 \param dir The param_common structure
 \param samples_struct The samples structure
 \param pos_param The param_sanePos structure
 \param proc_param The param_sanePre structure
 \param structPS The param_sanePS structure
 \param sanePic_struct The param_sanePic structure
 \param saneInv_struct The param_saneInv structure
 \param fcut_file A string containing the fcut filename
 \param check_struct The saneCheck structure
 \param size Total number of processors, in case paraframe or parabolo is defined
 \param rank The processor rank given by MPI_Comm_rank, in case paraframe or parabolo is defined
 \return A flag corresponding to an error code, or 0
 */
uint16_t parse_saneCheck_ini_file(char * ini_name, string &output, struct param_common &dir,
		struct samples &samples_struct, struct param_sanePos &pos_param, struct param_sanePre &proc_param,
		struct param_sanePS &structPS, struct param_saneInv &saneInv_struct, struct param_sanePic &sanePic_struct, struct saneCheck &check_struct, int rank, int size);

//! Reads saneCheck structure informations in ini file
/*!
 \param check_struct The saneCheck structure
 \param ini A pointer to a dictionnary opened using iniparser library
 */
void read_saneCheck_ini(dictionary	*ini , struct saneCheck &check_struct);

//! Print saneCheck structure informations to screen
/*!
 \param check_struct The saneCheck structure
 */
void print_saneCheck_ini(struct saneCheck check_struct);

#define PARSE_SANECHECK_H_


#endif /* PARSE_SANECHECK_H_ */
