#include "struct_definition.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#ifndef PARSE_SANECHECK_H_

uint16_t parse_saneCheck_ini_file(char * ini_name, string &output, struct param_common &dir,
		struct samples &samples_struct, struct param_sanePos &pos_param, struct param_sanePre &proc_param,
		struct param_sanePS &structPS, struct param_saneInv &saneInv_struct, struct param_sanePic &sanePic_struct, struct saneCheck &check_struct, int rank, int size);

int read_saneCheck_ini(dictionary	*ini , struct saneCheck &check_struct);
void print_saneCheck_ini(struct saneCheck check_struct);

#define PARSE_SANECHECK_H_


#endif /* PARSE_SANECHECK_H_ */
