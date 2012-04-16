#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include "utilities.h"

using namespace std;




char* tableFormat(std::vector<string> strings) {
	// Return the table Format for the given vector of string
	stringstream stream_tform;
	string string_tform;
	char* p_tform;

	stream_tform << maxStringLength(strings) << "A";
	string_tform = stream_tform.str();

	p_tform = new char[string_tform.size() + 1];
	strcpy(p_tform, string_tform.c_str());
	return p_tform;
}

long maxStringLength(std::vector<string> strings) {
	// Return the longest string length of a string vector
	unsigned long maxSize = 0;
	std::vector<string>::iterator itString;

	for (itString = strings.begin(); itString != strings.end(); itString++) {
		string iString = *(itString);
		if (iString.size() > maxSize)
			maxSize = iString.size();

	}
	return maxSize;
}

char** vString2carray(std::vector<string> strings) {
	// Transform a vector of string into a array of char

	int stringLength = maxStringLength(strings);
	int nBolos = strings.size();

	char **data;

	data = new char*[nBolos];

	for (int i = 0; i < nBolos; i++) {
		data[i] = new char[stringLength+1]; // TODO : +\0 ?
		strcpy(data[i], strings[i].c_str());
	}

	return data;

}

void vDouble2carray(std::vector<double> doubles, double ** data, long *size){
	// Transform a vector of double into a array of double, safe way

	*size = doubles.size();
	*data = new double[(*size)];
	for (long ii=0; ii<(*size); ii++){
		(*data)[ii] = doubles[ii];
	}

}

