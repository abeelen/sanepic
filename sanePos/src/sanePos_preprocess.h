/*
 * sanePre_preprocess.h
 *
 *  Created on: 25 juin 2009
 *      Author: matthieu
 */

/*! \file restypedef.cpp
 *  sanePre_preprocess.h
 *  Created on: 25 juin 2009
 *      Author: matthieu
 */

#ifndef SANEPOS_PREPROCESS_H_
#define SANEPOS_PREPROCESS_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <cstdlib>
#include <cmath>

#include <fcntl.h>
#include <unistd.h>

#include <vector>
#include <stdio.h>
#include <string>
#include <algorithm>


#include "sanePos_map_making.h"
#include "boloIO.h"
#include "dataIO.h"
#include "blastSpecific.h"
#include "inline_IO2.h"

using namespace std;

/*!
 *  find_coordinates_in_map : Output : ra_min, ra_max, dec_min, dec_max
 *
 */
void find_coordinates_in_map(long ndet,std::vector<string> bolonames,string bextension,
		string fextension,string file_offsets,foffset *foffsets,float *scoffsets,	/*double *offsets,*/long iframe_min, long iframe_max,
		long *fframes,long *nsamples,string dirfile,string ra_field,string dec_field,string phi_field, string scerr_field,
		string flpoint_field,int nfoff,double pixdeg, int *&xx, int *&yy,int nn, double *&coordscorner, double *tancoord,
		double *tanpix, bool bfixc, double radius, double *offmap, double *srccoord, char type, double *&ra,double *&dec,double *&phi,double *&scerr, unsigned char *&flpoint,double &ra_min,double &ra_max,double &dec_min,double &dec_max);

/*!
 *  compute_indpsrc_addnpix
 *  Input : nn : size of the map in pixels
 * ntotscan : total number of scans
 *  xxi, xxf, yyi, yyf : box for crossing constraint removal corner coordinates
 *  indpsrc : pixel indice of the pixel included in box for CCR
 *  npixsrc : total number of pix in box for CCR
 *  mask : 0 if the pixel is contained in box for CCR, 1 otherwise
 */
long Compute_indpsrc_addnpix(int nn, long ntotscan,std::vector<long> xxi, std::vector<long> xxf,
		std::vector<long> yyi, std::vector<long> yyf, long* &indpsrc, long &npixsrc, unsigned char* &mask);

/*!
 *  Get coordinates of pixels that are seen
 *  Compute the position to pixel projetcion matrices :
 *  One binary file per bolometer and per scan
 */
void compute_seen_pixels_coordinates(int ndet, long ntotscan,string outdir,std::vector<string> bolonames,string bextension, string fextension,string termin,
		string file_offsets,foffset *foffsets,float *scoffsets,long iframe_min, long iframe_max,long *fframes,
		long *nsamples,string dirfile,string ra_field,string dec_field,string phi_field, string scerr_field,
		string flpoint_field,int nfoff,double pixdeg, int *&xx, int *&yy,unsigned char *&mask, int &nn, double *&coordscorner, double *tancoord,
		double *tanpix, bool bfixc, double radius, double *offmap, double *srccoord, char type, double *&ra,double *&dec,
		double *&phi,double *&scerr, unsigned char *&flpoint,int shift_data_to_point,double &ra_min,double &ra_max,double &dec_min,double &dec_max,unsigned char* &flag,
		long napod, double errarcsec, bool NOFILLGAP,bool flgdupl,int factdupl, long addnpix, unsigned char *&rejectsamp, long *&samptopix, long *&pixon, int rank,
		long *indpsrc, long npixsrc, int &flagon,bool &pixout);

/*!
 *   compute indpix : pixels indices used for data projection/deprojection.
 *   npix is the number of filled pixels
 *   Pixon : indicates if the pixel is projected on the map or not
 *   addnpix : number of pixels to add : depends on box for CCR and duplication factor
 */
int compute_indpix(long *&indpix,int factdupl,int nn,long addnpix, long* pixon);

int map_offsets(string file_frame_offsets, long ntotscan, float *&scoffsets, foffset *&foffsets, long *fframes, int rank);

#endif /* SANEPOS_PREPROCESS_H_ */
