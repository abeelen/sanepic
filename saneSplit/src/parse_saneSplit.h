#ifndef PARSE_SANESPLIT_H_
#define PARSE_SANESPLIT_H_

#include "struct_definition.h"

//! Parse user command line and ini file
/*!
 * Only reads dir structure in ini file (part common)
 \param ini_name ini file name
 \param output The parser error string
 \param dir Output directory pathname
 \return An integer >0 in case ini file was not found, 0 otherwise
 */
int parse_saneSplit_ini_file(char * ini_name, std::string &output, struct param_common &dir);

#endif /* PARSE_SANESPLIT_H_ */
