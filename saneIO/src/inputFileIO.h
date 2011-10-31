#ifndef INPUTFILEIO_H_
#define INPUTFILEIO_H_

#include <string>
#include <vector>
#include <fstream>

#include "struct_definition.h"

//! Reads a detector list in an ascii file
/*!
 \param fname The file which contains the channel list
 \param bolos A vector in which the channel list is stored
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_strings(std::string fname, std::vector<std::string>& bolos);

//! Reads a double list in an ascii file
/*!
 \param fname The file which contains the channel list
 \param array An array in which the double list is stored
 \param size "array" size
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_double(std::string fname, std::vector<double> & array );

//! Given a string str, replace a part of this string (tobe_replace) with another (with_this)
/*!
 \param str A string that needs to be modified
 \param tobe_replace The part of "str" that has to be replace
 \param with_this The string that gonna replace "tobe_replace" in "str"
 \return The modified string
 */
std::string replace_all(std::string str, std::string tobe_replace, std::string with_this);

//! Removes from the input string both path and ".fits" extension
/*!
 \param path A string containing an input fits file filename
 \return The modified string
 */
std::string Basename(std::string path);
std::string FitsBasename(std::string path);

//! Removes from the input string ".fits" extension and replace all '.' with '_' character
/*!
 \param path A string containing an input fits file filename
 \return The modified string
 */
std::string dirfile_Basename(std::string path);

//! Get the number of samples for each scan contained in "inputFiles" vector
/*!
 \param inputdir Input directory path
 \param inputFiles A vector of string containing input fits filenames
 \param nsamples An array in which the number of samples for each scan is stored
 */
void readFrames(std::vector<std::string> &inputFiles, std::vector<long> &nsamples);

//! Given an input fits "filename", get its number of samples (called by readFrames routine)
/*!
 \param filename A string containing an input fits file global filename
 \return The number of samples
 */
long readFitsLength(std::string filename);

//! Given an input bololist filename (fname), Reads the channel list
/*!
 \param fname An input bololist filename
 \param bolonames This list is stored in this vector of string
 \param output The parser error string
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_channel_list(std::string &output, std::string fname, std::vector<std::string> &bolonames);

//! Given an input fitsfile list (fname), Reads the scan list (and processor order, if any)
/*!
 \param fname The input fitsfile list
 \param samples_struct A samples structure in which this file informations are stored
 \param output The parser error string
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
uint16_t read_fits_list(std::string &output, std::string fname, struct samples &samples_struct);

void skip_comment(std::ifstream &file, std::string &line);

#endif /* INPUTFILEIO_H_ */
