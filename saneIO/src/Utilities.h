/*
 * utilities.h
 *
 *  Created on: 24 janv. 2011
 *      Author: matthieu
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <utility>

#include <limits>

#include <stdint.h>

using namespace std;

#ifdef USE_MPI
#include "mpi.h"
#endif

template<typename T, size_t N> T * end(T (&ra)[N]) {
	return ra + N;
}

#define Elements_in(arrayname) (sizeof arrayname/sizeof *arrayname)


//! Given a variable "object", which type can be boolean, numerical or string, transform it as a string using stringstream
/*!
 \param object An object that can be a string (or char*), a boolean, or a numerical value
 \return The same content but casted as a string
 */
template <typename T>
std::string StringOf(const T& object){
	std::ostringstream os;
	os.precision(std::numeric_limits< double >::digits10);
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


/** Sorts a vector and returns index of the sorted values
 * \param Index Contains the index of sorted values in the original vector
 * \param data The vector to be sorted
 */
template<class T> void paired_sort(const std::vector<T> & data, std::vector<size_t> & Index ){
	// A vector of a pair which will contain the sorted value and its index in the original array
	std::vector<pair<T,size_t> > IndexedPair;
	IndexedPair.resize(data.size());
	for(size_t i=0;i<IndexedPair.size();++i)
	{
		IndexedPair[i].first = data[i];
		IndexedPair[i].second = i;
	}
	sort(IndexedPair.begin(),IndexedPair.end());
	Index.resize(data.size());
	for(size_t i = 0; i < Index.size(); ++i) Index[i] = IndexedPair[i].second;
}

/** Reorder a vector given an index vector
 *  From http://stackoverflow.com/questions/838384/reorder-vector-using-a-vector-of-indices
 *  \param data the target vector
 *  \order the ordered index vector
 */
template <typename T> void reorder( std::vector<T> & data, std::vector<std::size_t> const & order )
{
	std::vector<T> tmp;         // create an empty vector
	tmp.reserve( data.size() ); // ensure memory and avoid moves in the vector
	for ( std::size_t i = 0; i < order.size(); ++i ) {
		tmp.push_back( data[order[i]] );
	}
	data.swap( tmp );          // swap vector contents
}


uint16_t fill_channel_list(std::string &output, struct samples &samples_struct, int rank, int size);

//! Computes the size of each string contained in "str_vect"
/*!
 \param str_vect A vector containing a list of channels (their names)
 \return size_buff : A long that corresponds to maximum string size
 */
unsigned long max_size_str_vect(std::vector<std::string> str_vect);

long getTotalSystemMemory();
long getAvailableSystemMemory();

// Helper function for screen print
std::string prettyPrintSize(double size);
std::string prettyPrintAngle(float angle);

#ifdef USE_MPI

uint16_t MPI_Bcast_vector_string(std::vector<std::string> &fitsvect, int rank, MPI_Comm Comm);
uint16_t MPI_Bcast_vector_int(std::vector<int> & intvect, int rank, MPI_Comm Comm);
uint16_t MPI_Bcast_vector_long(std::vector<long> & longvect, int rank, MPI_Comm Comm);
uint16_t MPI_Bcast_dmatrix(double ** dmatrix, long nrow, long ncol, int root, MPI_Comm Comm);
uint16_t MPI_Reduce_dmatrix(   double **sendrecvbuf, long nrow, long ncol, MPI_Op op, int root, MPI_Comm Comm);
uint16_t MPI_Allreduce_dmatrix(double **sendrecvbuf, long nrow, long ncol, MPI_Op op, MPI_Comm Comm);

#endif


#endif /* UTILITIES_H_ */
