/*
 * inline_IO2.h
 *
 *  Created on: 23 juin 2009
 *      Author: matthieu
 */

#ifndef INLINE_IO2_H_
#define INLINE_IO2_H_


#include <fftw3.h>


extern "C" {
#include "nrutil.h"
#include "getdata.h"
}


// sanePos functions

bool compute_dirfile_format_file(std::string outdir, struct detectors det, long ntotscan, int rank);
bool compute_dirfile_format_noisePS(std::string outdir, std::vector<std::string> det, std::string suffix);
bool compute_dirfile_format_fdata(std::string outdir, struct detectors det, long ntotscan, int rank);

/*!
 * Writes indpix in a binary file \n
 * -ind_size is the indpix readed size \n
 * -npix the number of filled pixels \n
 * -indpix the pixels indices in the map \n
 * -flagon : boolean, indicates whether a sample has been rejected or not
 */
void write_indpix(long long ind_size, long long npix, long long *indpix,  std::string outdir, int flagon);
void read_indpix(long long &ind_size, long long &npix, long long *&indpix, std::string outdir, int &flagon);

/*!
 * Reads indpix from a binary file \n
 * -ind_size is the indpix readed size \n
 * -npix the number of filled pixels \n
 * -indpix the pixels indices in the map \n
 * -flagon : boolean, indicates whether a sample has been rejected or not
 */

void write_indpsrc(long long map_size,  long long  npixsrc, long long * indpsrc,   std::string outdir);
void  read_indpsrc(long long &map_size, long long &npixsrc, long long *&indpsrc,   std::string outdir);



/*!
 * Writes PNd in a binary file \n
 * -PNd is the Projected Noised Data \n
 * -npix the number of filled pixels \n
 */
void write_PNd(double *PNd, long long npix, std::string outdir);

/*!
 * Reads PNd from a binary file \n
 * -PNd is the Projected Noised Data \n
 * -npix the number of filled pixels \n
 */
void read_PNd(double *&PNdtot, long long &npix, std::string outdir);

/*!
 * Writes samptopix in a binary file \n
 * -ns is the number of samples for the iframe frame \n
 * -iframe is the frame number \n
 * -idet is the detector number \n
 * -samptopix is sample to pixel projection matrix
 */
void write_samptopix(long ns, long long *&samptopix,  std::string outdir, long iframe, std::string boloname);

/*!
 * Reads samptopix from a binary file \n
 * -ns is the number of samples for the iframe frame \n
 * -iframe is the frame number \n
 * -idet is the detector number \n
 * -samptopix is sample to pixel projection matrix
 */
void read_samptopix(long ns, long long *&samptopix, std::string outdir, long iframe, std::string boloname);

/*!
 * Writes fourier transform of the data in a binary file \n
 * -ns is the number of samples for the iframe frame \n
 * -iframe is the frame number \n
 * -idet is the detector number \n
 * -fdata is the fourier transform of the iframe data, for detector idet
 */
void write_fdata(long ns, fftw_complex *fdata, std::string outdir, long idet, long iframe, std::vector<std::string> bolonames);

/*!
 * Reads PS for one bolo in a binary file \n
 * -nbins is the number of bins \n
 * -ell is the bins values \n
 * -Spn_all is the detector noise PS \n
 * -nameSpfile is the name of the binary file
 * -ndet = number of detector
 */
//void read_noise_file(long &nbins, double *&ell, double **&SpN_all, string nameSpfile, long ndet);

/*!
 * Reads fourier transform of the data from a binary file \n
 * -ns is the number of samples for the iframe frame \n
 * -iframe is the frame number \n
 * -idet is the detector number \n
 * -fdata is the fourier transform of the iframe data, for detector idet
 */
void read_fdata(long ns, fftw_complex *&fdata, std::string prefixe, std::string outdir, long idet, long iframe, std::vector<std::string> bolonames);

/*!
 * Writes fourier transform of the data in a binary file \n
 * -ns is the number of samples for the iframe frame \n
 * -iframe is the frame number \n
 * -idet is the detector number \n
 * -fdata is the fourier transform of the iframe data, for detector idet
 */
void write_fPs(long ns, fftw_complex *fdata, std::string outdir, long idet, long iframe, std::vector<std::string> bolonames);

/*!
 * Writes information for a second run of sanepic \n
 * -nn is the size of a map side
 * -npix is the number of filled pixels
 * -pixdeg is the pixel size in degrees
 * -tancoord is the tangent point coordinates
 * -tanpix are the tangent point pixel indices
 * -coordsyst : the coordinate system
 * -flagon indicates whether a sample was rejected
 * -indpix is the pixel indices in the map
 */
//void write_info_for_second_part(string outdir, int NAXIS1, int NAXIS2, int npix,
//		double pixdeg, double *tancoord, double* tanpix, int coordsyst, bool flagon, long* indpix);

/*!
 * Reads mixing matrix in a .txt file \n
 * -MixMatfile is the name of the mixing matrix file
 * -ndet number of detector
 * -ncomp2 = number of component
 * -mixmat is the mixing matrix
 */
void read_mixmat_txt(std::string MixMatfile, long ndet, long ncomp, double **&mixmat);



#endif
