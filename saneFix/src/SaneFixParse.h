/*
 * SaneFixParse.h
 *
 */

#ifndef SANEFIXPARSE_H_
#define SANEFIXPARSE_H_

#include "StructDefinition.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

uint32_t parse_saneFix_ini_file(char * ini_name, string &output, struct param_common &dir, 		struct samples &samples_struct, struct param_sanePos &pos_param, struct param_saneProc &proc_param,
		struct param_sanePS &structPS, struct param_saneInv &saneInv_struct, struct param_sanePic &sanePic_struct, struct param_saneFix &Fix_struct, int size, int rank);




#endif /* SANEFIXPARSE_H_ */
