/*
 * Corr_preprocess.h
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */

#ifndef CORR_PREPROCESS_H_
#define CORR_PREPROCESS_H_


#include "todprocess.h"
#include "map_making.h"



void do_PtNd(double *PNd, string *extentnoiseSp_all, string noiseSppreffile,
		string dir, string prefixe, string termin, string *bolonames,
		double f_lppix, double fsamp, long ff, long ns, long ndet, int size,
		int rank, long *indpix, long nn, long npix, long iframe, double *Mp, long *hits);




void write_ftrProcesdata(double *S, long *indpix, long *indpsrc, int nn, long npix,
		long npixsrc, long ntotscan, long addnpix, bool flgdupl, int factdupl,
		int fillg, string dir, string termin, double errarcsec, string dirfile,
		string scerr_field, string flpoint_field, string *bolonames,
		string bextension, string fextension, /*string cextension,*/
		int shift_data_to_point, double f_lppix, long ff, long ns,
		long napod, long ndet, bool NORMLIN, bool NOFILLGAP,
		long iframe);



void write_tfAS(double *S, long *indpix, int nn, long npix, bool flgdupl, int factdupl,
		string dir, string termin, long ff, long ns, long ndet, long iframe);

#endif /* CORR_PREPROCESS_H_ */