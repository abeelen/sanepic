#include "StructDefinition.h"

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
 \param proc_param The param_saneProc structure
 \param structPS The param_sanePS structure
 \param sanePic_struct The param_sanePic structure
 \param saneInv_struct The param_saneInv structure
 \param check_struct The saneCheck structure
 \param size Total number of processors, in case paraframe or parabolo is defined
 \param rank The processor rank given by MPI_Comm_rank, in case paraframe or parabolo is defined
 \return A flag corresponding to an error code, or 0
 */
uint32_t parse_saneCheck_ini_file(char * ini_name, string &output, struct param_common &dir,
		struct samples &samples_struct, struct param_sanePos &pos_param, struct param_saneProc &proc_param,
		struct param_sanePS &structPS, struct param_saneInv &saneInv_struct, struct param_sanePic &sanePic_struct, struct param_saneCheck &check_struct, int size, int rank);

#define PARSE_SANECHECK_H_


#endif /* PARSE_SANECHECK_H_ */
