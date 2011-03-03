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

template <typename T>
std::string StringOf(const T& object){
	std::ostringstream os;
	os << object;
	return os.str();
}


/*! Returns the table Format for the given vector of string*/
char* tableFormat(std::vector<std::string> strings);

/*! Return the longest string length of a string vector */
long maxStringLength(std::vector<std::string> strings);

/*! Transform a vector of string into a array of char*/
char** vString2carray(std::vector<std::string> strings);


#endif /* UTILITIES_H_ */
