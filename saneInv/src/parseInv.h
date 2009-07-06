/*
 * parseInv.h
 *
 *  Created on: 18 juin 2009
 *      Author: matthieu
 */

#ifndef PARSEINV_H_
#define PARSEINV_H_
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <list>
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
//int parse_sanePos_ini_file(char * ini_name);
int parse_saneInv_ini_file(char * ini_name, string &fname,string &boloname, string &noiseSp_dir_output, string &extentnoiseSp);

#endif /* PARSEINV_H_ */
