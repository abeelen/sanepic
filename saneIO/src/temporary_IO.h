#ifndef INLINE_IO2_H_
#define INLINE_IO2_H_


#include <fftw3.h>
#include "struct_definition.h"


extern "C" {
#include "nrutil.h"
}

// dirfile functions
int write_data_flag_to_dirfile(struct param_common dir, struct samples samples_struct, long iframe_min, long iframe_max);
int read_data_from_dirfile(std::string tmp_dir, std::string filename, std::string field, double *&data);
int read_flag_from_dirfile(std::string tmp_dir, std::string filename, std::string field, int *&mask);
// dirfile functions
//bool compute_dirfile_format_file(std::string tmp_dir, struct dirfile_fragment &dirf);
//bool compute_dirfile_format_noisePS(std::string outdir, std::vector<std::string> det, std::string suffix);
//bool compute_dirfile_format_fdata(std::string outdir, struct samples samples_struct, int rank);
//bool compute_dirfile_format_fdata(std::string outdir, struct samples samples_struct, int rank);

/*!
 * Writes indpix in a binary file \n
 * -ind_size is the indpix readed size \n
 * -npix the number of filled pixels \n
 * -indpix the pixels indices in the map \n
 * -flagon : boolean, indicates whether a sample has been rejected or not
 */
int write_indpix(long long ind_size, long long npix, long long *indpix,  std::string outdir, int flagon);
int read_indpix(long long &ind_size, long long &npix, long long *&indpix, std::string outdir, int &flagon);

/*!
 * Reads indpix from a binary file \n
 * -ind_size is the indpix readed size \n
 * -npix the number of filled pixels \n
 * -indpix the pixels indices in the map \n
 * -flagon : boolean, indicates whether a sample has been rejected or not
 */

int write_indpsrc(long long map_size,  long long  npixsrc, long long * indpsrc,   std::string outdir);
int  read_indpsrc(long long &map_size, long long &npixsrc, long long *&indpsrc,   std::string outdir);



/*!
 * Writes PNd in a binary file \n
 * -PNd is the Projected Noised Data \n
 * -npix the number of filled pixels \n
 */
int write_PNd(double *PNd, long long npix, std::string outdir);

/*!
 * Reads PNd from a binary file \n
 * -PNd is the Projected Noised Data \n
 * -npix the number of filled pixels \n
 */
int read_PNd(double *&PNdtot, long long &npix, std::string outdir);

/*!
 * Writes samptopix in a binary file \n
 * -ns is the number of samples for the iframe frame \n
 * -iframe is the frame number \n
 * -idet is the detector number \n
 * -samptopix is sample to pixel projection matrix
 */
int write_samptopix(long ns, long long *&samptopix,  std::string tmpdir, std::string filename, std::string boloname);

/*!
 * Reads samptopix from a binary file \n
 * -ns is the number of samples for the iframe frame \n
 * -iframe is the frame number \n
 * -idet is the detector number \n
 * -samptopix is sample to pixel projection matrix
 */
int read_samptopix(long ns, long long *&samptopix, std::string tmpdir, std::string filename, std::string boloname);

/*!
 * Writes fourier transform of the data in a binary file \n
 * -ns is the number of samples for the iframe frame \n
 * -iframe is the frame number \n
 * -idet is the detector number \n
 * -fdata is the fourier transform of the iframe data, for detector idet
 */
int write_fdata(long ns, fftw_complex *fdata, std::string prefixe, std::string tmpdir, long idet, std::string filename, std::vector<std::string> bolonames);

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
int read_fdata(long ns, fftw_complex *&fdata, std::string prefixe, std::string tmpdir, long idet, std::string filename, std::vector<std::string> bolonames);

/*!
 * Reads mixing matrix in a .txt file \n
 * -MixMatfile is the name of the mixing matrix file
 * -ndet number of detector
 * -ncomp2 = number of component
 * -mixmat is the mixing matrix
 */
int read_mixmat_txt(std::string MixMatfile, long ndet, long ncomp, double **&mixmat);



#endif
