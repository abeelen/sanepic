/*
 * psdIO.h
 *
 *  Created on: 16 juin 2009
 *      Author: matthieu
 */

#ifndef PSDIO_H_
#define PSDIO_H_


#include <iostream>
#include <string>

extern "C" {
#include <fitsio.h>
}

using namespace std;

int write_psd_tofits(string fname, long nx, long ny,char dtype, void * psd1d);
#endif /* PSDIO_H_ */
