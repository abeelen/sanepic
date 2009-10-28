/*
 * sanepic_preprocess.h
 *
 *  Created on: 2 juil. 2009
 *      Author: matthieu
 */

#ifndef SANEPIC_PREPROCESS_H_
#define SANEPIC_PREPROCESS_H_

#include <vector>

using namespace std;

void sanepic_preprocess(int NAXIS1, int NAXIS2, std::vector<long> xxi, std::vector<long> xxf,
		std::vector<long> yyi, std::vector<long> yyf, long *&indpsrc, long &npixsrc,
		long ntotscan, long &addnpix);

#endif /* SANEPIC_PREPROCESS_H_ */
