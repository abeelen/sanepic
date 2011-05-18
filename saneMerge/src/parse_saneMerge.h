#ifndef PARSE_SANEMERGE_H_
#define PARSE_SANEMERGE_H_

//! Parse user command line
/*!
 \param opt_name User command line : argv
 \param output The parser error string
 \param arg Number of argument in user command line : argc
 \param dir Output directory pathname
 \param samples_struct The samples structure
 \return An integer corresponding to an error code, or 0 if everything went OK
 */
int parse_saneMerge_ini_file(char * opt_name[], std::string &output, int arg, struct param_common &dir,
		struct samples &samples_struct);


#endif /* PARSE_SANEMERGE_H_ */
