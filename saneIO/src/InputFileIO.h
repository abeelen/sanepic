#ifndef INPUTFILEIO_H_
#define INPUTFILEIO_H_

#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include <stdint.h>

#include "ErrorCode.h"
#include "StructDefinition.h"

// word_count from http://www.cplusplus.com/forum/general/30929/
// general case, stream interface
inline size_t word_count(std::istream& is)  // can pass an open std::ifstream() to this if required
{
	size_t c = 0;
	for(std::string w; is >> w; ++c);
	return c;
}
// simple string interface
inline size_t word_count(const std::string& str)
{
	std::istringstream iss(str);
	return word_count(iss);
}

// Helper for ASCII file IO
// Template function needs to be declared once as they are generated at compilation...
uint16_t read_file_line(std::string &output, std::string fname, std::vector<std::string> & content );

template <typename T> uint16_t read_file(std::string & output, std::string fname, std::vector<T> & output_T ) {
  /*
   * Read an ASCII file and return one column
   */

  std::vector<std::string> file_content;

  T i_T;

  if ( read_file_line(output, fname, file_content) )
    return FILE_PROBLEM;

  output_T.clear();

  // Fill output_string & output_float
  for (std::vector<std::string>::iterator it = file_content.begin(); it!=file_content.end(); ++it) {
    // Check for a least two columns
    if ( word_count(*it) < 1 ){
      output += "EE - "+fname+" must contain at least one column on each line\n";
      return FILE_PROBLEM;
    }
    // stream the columns...
    std::istringstream iline(*it);
    iline >> i_T;
    output_T.push_back(i_T);
  }

  return 0;

}

template <typename T, typename U> uint16_t read_file_2col(std::string &output, std::string fname, std::vector<T> & output_T, std::vector<U> & output_U ) {
  /**
   * Read an ASCII file and return two columns depending on the desired type
   */

  std::vector<std::string> file_content;

  T i_T;
  U i_U;

  if ( read_file_line(output, fname, file_content) )
    return FILE_PROBLEM;

  output_T.clear();
  output_U.clear();

  // Fill output_string & output_float
  for (std::vector<std::string>::iterator it = file_content.begin(); it!=file_content.end(); ++it) {
    // Check for a least two columns
    if ( word_count(*it) < 2 ){
      output += "EE - "+fname+" must contain at least 2 column on each line\n";
      return FILE_PROBLEM;
    }
    // stream the columns...
    std::istringstream iline(*it);
    iline >> i_T >> i_U;
    output_T.push_back(i_T);
    output_U.push_back(i_U);
  }

  return 0;

}


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
void readFramesFromFits(struct samples &samples_struct);

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
int readChannelList(std::string &output, std::string fname, std::vector<std::string> &bolonames);

//! Given an input fitsfile list (fname), Reads the scan list (and processor order, if any)
/*!
 \param fname The input fitsfile list
 \param samples_struct A samples structure in which this file informations are stored
 \param output The parser error string
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
uint16_t readFitsList(std::string &output, std::string fname, struct samples &samples_struct);


#endif /* INPUTFILEIO_H_ */
