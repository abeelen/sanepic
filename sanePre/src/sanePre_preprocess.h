/*
 * sanePre_preprocess.h
 *
 *  Created on: 25 juin 2009
 *      Author: matthieu
 */

#ifndef SANEPRE_PREPROCESS_H_
#define SANEPRE_PREPROCESS_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <cstdlib>

#include <fcntl.h>
#include <unistd.h>
//#include <list>
#include <vector>
#include <stdio.h>
#include <string>
#include <algorithm>

#include "todprocess.h"
#include "map_making.h"


long Compute_indpsrc_addnpix(int nn, long ntotscan,std::vector<long> xxi, std::vector<long> xxf,
		std::vector<long> yyi, std::vector<long> yyf, long* &indpsrc, long &npixsrc);


#endif /* SANEPRE_PREPROCESS_H_ */
