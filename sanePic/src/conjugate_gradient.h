/*
 * conjugate graident.h
 *
 *  Created on: 3 juil. 2009
 *      Author: matthieu
 */

#ifndef CONJUGATEGRAIDENT_H_
#define CONJUGATEGRAIDENT_H_


#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <sstream>
#include "todprocess.h"
#include "map_making.h"

//#include "sane_io.h"
#include "binaryFileIO.h"
#include "boloIO.h"
#include "dataIO.h"
#include "imageIO.h"
#include "inline_IO2.h"
//#include "mpi_architecture_builder.h"

#include "parseSanepic.h"
#include "sanepic_preprocess.h"

//#include "estimPS_sanepic.h"
#include "Corr_preprocess.h"
#include "NoCorr_preprocess.h"
#include "mpi_architecture_builder.h"
#include <time.h>
#include <fftw3.h>
//#include <fcntl.h>
//#include <unistd.h>
#include <list>
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;

void sanepic_conjugate_gradient(bool flgdupl, int npix, double* &S, long iframe_min, long iframe_max,
		long *nsamples, long *fframes,std::vector<double> fcut,double f_lp,double fsamp,
		long *indpix, int nn, int factdupl, string poutdir, string termin, string termin_internal, long ndet,
		string *extentnoiseSp_all,string noiseSppreffile, std::vector<string> bolonames, int size_det,
		int rank_det, int iterw, double pixdeg,double *tancoord, double *tanpix,int coordsyst,
		long *indpsrc, long npixsrc, int flagon, bool projgaps, int rank, bool CORRon,
		string dirfile, double *&PNdtot, long ntotscan,long addnpix,bool NORMLIN,bool NOFILLGAP,
		long napod,int shift_data_to_point,bool remove_polynomia,string fextension,string bextension,
		string flpoint_field,string scerr_field, string outdir);

#endif /* CONJUGATEGRAIDENT_H_ */
