/*
 * parse_FBFO.h
 *
 *  Created on: 15 juin 2009
 *      Author: matthieu
 */

#ifndef PARSE_FBFO_H_
#define PARSE_FBFO_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include "boloIO.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

using namespace std;

int parse_FBFO(char * ini_name, string &fname, long &ntotscan, std::vector<long> &fframes_vec, std::vector<long> &nsamples_vec);

#endif /* PARSE_FBFO_H_ */
