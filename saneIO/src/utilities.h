/*
 * utilities.h
 *
 *  Created on: 24 janv. 2011
 *      Author: matthieu
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>


//! Given a variable "object", which type can be boolean, numerical or string, transform it as a string using stringstream
/*!
 \param object An object that can be a string (or char*), a boolean, or a numerical value
 \return The same content but casted as a string
 */
template <typename T>
std::string StringOf(const T& object){
	std::ostringstream os;
	os << object;
	return os.str();
}


//! Returns the cfitsio table Format for a given vector of string
/*!
 \param strings A vector of strings
 \return A char* containing the tform
 */
char* tableFormat(std::vector<std::string> strings);

//! Return the longest string length of a string vector
/*!
 \param strings A vector of strings
 \return The longest string length as a long
 */
long maxStringLength(std::vector<std::string> strings);

//! Transform a vector of string into a array of char
/*!
 \param strings A vector of strings
 \return A char** containing the same informations as strings, but using a char* array format
 */
char** vString2carray(std::vector<std::string> strings);


//! Transform a vector of double into a array of double
/*!
 \param doubles A vector of double
 \return A double * array and long containing the same informations as strings, but using a char* array format
 */
void vDouble2carray(std::vector<double> doubles, double ** data, long *size);

#endif /* UTILITIES_H_ */
