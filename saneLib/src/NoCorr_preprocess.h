/*
 * NoCorr_preprocess.h
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */

#ifndef NOCORR_PREPROCESS_H_
#define NOCORR_PREPROCESS_H_

#include <vector>
#include <string>

#include "todprocess.h"
#include "map_making.h"
#include "inline_IO2.h"

void do_PtNd_nocorr(double *PNd, string *extentnoiseSp_all, string noiseSppreffile,
		string dir, /* string termin, double errarcsec,*/ string dirfile,
		/*string scerr_field, string flpoint_field,*/ std::vector<string> bolonames, string *fits_table,
		/*string bextension, string fextension,*//* string cextension,*/
		int shift_data_to_point, double f_lppix, double f_lppix_Nk,
		double fsamp, long ntotscan, long addnpix, bool flgdupl, int factdupl,
		int fillg, long ff, long ns, long napod, long ndet,
		long *indpix, long *indpsrc, long nn, long npix, long npixsrc,
		bool NORMLIN, bool NOFILLGAP,bool remove_polynomia,  long iframe, double *S);



void do_PtNPS_nocorr(double *S, string *extentnoiseSp_all, string noiseSppreffile, string dir,
		/*string termin,*/ string dirfile, std::vector<string> bolonames, double f_lppix,
		double fsamp, bool flgdupl, int factdupl, long ff, long ns,
		long ndet, long *indpix, long nn, long npix,
		long iframe, double *PtNPmatS, double *Mp, long *hits);

#endif /* NOCORR_PREPROCESS_H_ */
