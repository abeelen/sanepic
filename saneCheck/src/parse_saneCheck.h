#include "struct_definition.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#ifndef PARSE_SANECHECK_H_

int parse_saneCheck_ini_file(char * ini_name, struct param_common &dir,
		std::vector<detectors> &detector_tab,struct samples &samples_struct, double &fsamp, struct saneCheck &check_struct, int rank);

int read_saneCheck_ini(dictionary	*ini , struct saneCheck &check_struct, int rank);
void print_saneCheck_ini(struct saneCheck check_struct, int rank);

#define PARSE_SANECHECK_H_


#endif /* PARSE_SANECHECK_H_ */
