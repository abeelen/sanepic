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

#include "mpi_architecture_builder.h"

using namespace std;

/*!
 *  find_coordinates_in_map : Output : ra_min, ra_max, dec_min, dec_max
 *
 */

void find_coordinates_in_map(long ndet,std::vector<string> bolonames, string *fits_table,/*string bextension,
		string fextension,*//*string file_offsets,foffset *foffsets,float *scoffsets,	double *offsets,*/long iframe_min, long iframe_max,
		/*long *fframes,*/long *nsamples,string dirfile,/*string ra_field,string dec_field,string phi_field, string scerr_field,
		string flpoint_field,int nfoff, */double pixdeg, int *&xx, int *&yy,int nn, double *&coordscorner, double *tancoord,
		double *tanpix, bool bfixc, double radius, /*double *offmap,*/ double *srccoord, char type, double *&ra,double *&dec,double *&phi, short *&flpoint,double &ra_min,double &ra_max,double &dec_min,double &dec_max,bool default_projection);

/*!
 *  compute_indpsrc_addnpix
 *  Input : nn : size of the map in pixels
 * ntotscan : total number of scans
 *  xxi, xxf, yyi, yyf : box for crossing constraint removal corner coordinates
 *  indpsrc : pixel indice of the pixel included in box for CCR
 *  npixsrc : total number of pix in box for CCR
 *  mask : 0 if the pixel is contained in box for CCR, 1 otherwise
 */

void Compute_indpsrc_addnpix(unsigned long NAXIS1, unsigned long NAXIS2, long ntotscan, std::vector<struct box> boxFile,
		long &addnpix, long* &indpsrc, long &npixsrc, unsigned short* &mask);

/*!
 *  Get coordinates of pixels that are seen
 *  Compute the position to pixel projetcion matrices :
 *  One binary file per bolometer and per scan
 */

void compute_seen_pixels_coordinates(long ntotscan,string outdir,std::vector<string> bolonames, string *fits_table,/*string bextension, string fextension, string termin, */
		/*string file_offsets,foffset *foffsets,float *scoffsets, */ long iframe_min, long iframe_max,/*long *fframes,*/
		long *nsamples,string dirfile/*,string ra_field,string dec_field,string phi_field, string scerr_field,
		string flpoint_field,int nfoff*/,double pixdeg, int *&xx, int *&yy,unsigned short *&mask, int &nn, double *&coordscorner, double *tancoord,
		double *tanpix, bool bfixc, double radius, /*double *offmap,*/ double *srccoord, char type, double *&ra,double *&dec,
		double *&phi, short *&flpoint,int shift_data_to_point,double &ra_min,double &ra_max,double &dec_min,double &dec_max,short* &flag,
		long napod, double errarcsec, bool NOFILLGAP,bool flgdupl,int factdupl, long addnpix, unsigned char *&rejectsamp, long *&samptopix, long *&pixon, int rank,
		long *indpsrc, long npixsrc, int &flagon,bool &pixout);

void computePixelIndex(long ntotscan,string outdir, std::vector<string> bolonames,
		string *fits_table, long iframe_min, long iframe_max, unsigned long *nsamples,
		struct wcsprm & wcs, long NAXIS1, long NAXIS2,
		unsigned short *&mask,
		long napod,  bool NOFILLGAP,bool flgdupl, int factdupl,
		long addnpix, long *&pixon, int rank,
		long *indpsrc, long npixsrc, int &flagon, bool &pixout);

void computePixelIndex_HIPE(long ntotscan,string outdir, std::vector<string> bolonames,
		string *fits_table, long iframe_min, long iframe_max, unsigned long *nsamples,
		struct wcsprm & wcs, long NAXIS1, long NAXIS2,
		unsigned short *&mask,
		long napod,  bool NOFILLGAP,bool flgdupl, int factdupl,
		long addnpix, long *&pixon, int rank,
		long *indpsrc, long npixsrc, int &flagon, bool &pixout);


// TODO: Remove this
//int map_offsets(string file_frame_offsets, long ntotscan, float *&scoffsets, foffset *&foffsets, long *fframes, int rank);

#endif /* SANEPOS_PREPROCESS_H_ */
