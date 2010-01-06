
#ifndef PARSEPOS_H_
#define PARSEPOS_H_

#include <vector>
#include "mpi_architecture_builder.h"

/*!
 * - Parse sanePos input options using the ini file given in the command line \n
 * - Verify the good usage of each option, deal with incorrect values, warn the user when a line is \n
 * missing in the ini file !
 * - Initialise sanePos variable
 */
int parse_sanePos_ini_file(char * ini_name,struct param_process &com, struct param_positions &pos_param, struct directories &dir,
		struct detectors &det,struct samples &samples_struct,
		 int rank);

#endif /* PARSEPOS_H_ */
