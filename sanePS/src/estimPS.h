/*
 * estimPS.h
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */

#ifndef ESTIMPS_H_
#define ESTIMPS_H_


#include <vector>
#include <algorithm>
#include <sstream>
#include <time.h>

#include "todprocess.h"
#include "map_making.h"
#include "inline_IO2.h"
#include "covMatrixIO.h"
#include "psdIO.h"

#include "estimPS_steps.h"

//#include "sane_io.h"

extern "C" {
#include <fitsio.h>
}

void EstimPowerSpectra(double fsamp, long ns, long ff, long ndet, int NAXIS1, int NAXIS2, long npix, long napod,
		long iframe, bool flgdupl, int factdupl, long *indpix,
		double *S, string MixMatfile, std::vector<string> bolonames, string dirfile, string bextension,
		string fextension, int shift_data_to_point, string dir,
		bool NORMLIN, bool NOFILLGAP, bool remove_polynomia, string noiseSppreffile,
		string extentnoiseSp, string outdirSpN);


#endif /* ESTIMPS_H_ */
