/*
 * sanepic_preprocess.h
 *
 *  Created on: 2 juil. 2009
 *      Author: matthieu
 */

#ifndef SANEPIC_PREPROCESS_H_
#define SANEPIC_PREPROCESS_H_



#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdlib>

//#include "sane_io.h"
#include "binaryFileIO.h"
#include "dataIO.h"
#include "imageIO.h"
#include "inline_IO2.h"
//#include "mpi_architecture_builder.h"

#include "parseSanepic.h"
#include "sanepic_preprocess.h"

//#include "estimPS.h"
#include "Corr_preprocess.h"
#include "NoCorr_preprocess.h"
//#include "mpi_architecture_builder.h"
#include <time.h>
#include <fftw3.h>
//#include <fcntl.h>
//#include <unistd.h>
#include <list>
#include <stdio.h>
#include <stdlib.h>


using namespace std;

void sanepic_preprocess(int NAXIS1, int NAXIS2, std::vector<long> xxi, std::vector<long> xxf,
		std::vector<long> yyi, std::vector<long> yyf, long *&indpsrc, long &npixsrc,
		long ntotscan, long &addnpix);

#endif /* SANEPIC_PREPROCESS_H_ */
