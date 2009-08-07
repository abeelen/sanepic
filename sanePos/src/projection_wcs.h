/*
 * projection_wcs.h
 *
 *  Created on: 7 août 2009
 *      Author: matthieu
 */

#ifndef PROJECTION_WCS_H_
#define PROJECTION_WCS_H_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>

#include <sph.h>
#include <wcs.h>
#include <wcshdr.h>
#include <wcsfix.h>

extern "C" {
#include <fitsio.h>
}

void parser(struct wcsprm *wcs,double *tanpix,double *tancoord, double pixdeg);

void fits_header_generation(std::string outdir, const char *fits_file,double pixdeg, bool default_projection, double *tanpix,double *tancoord);

#endif /* PROJECTION_WCS_H_ */
