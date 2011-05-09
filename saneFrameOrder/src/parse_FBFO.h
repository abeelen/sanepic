#ifndef PARSE_FBFO_H_
#define PARSE_FBFO_H_


//! saneFrameOrder parser function
/*!
 * saneFrameOrder does not call parser_function from saneIO because it doesn't need all structure informations
 \param ini_name the ini file name given by user in the command line (usually argv[1])
 \param parser_output The parser error string
 \param dir The param_common structure
 \param samples_struct The samples structure
 \return An integer =-1 if there were an error, or 0 if everything went OK
 */
int parse_FBFO(char * ini_name, std::string &parser_output, struct samples &samples_struct,struct param_common &dir);

#endif /* PARSE_FBFO_H_ */
