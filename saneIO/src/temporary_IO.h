#ifndef INLINE_IO2_H_
#define INLINE_IO2_H_


#include <fftw3.h>
#include "struct_definition.h"


extern "C" {
#include "nrutil.h"
#include "getdata.h"
}

//! Reads data and flag tables in input scans and write them down to disk using getdata and dirfile format
/*!
  \param dir The param_common structure
  \param samples_struct The samples structure
  \param iframe_min Actual rank first frame indice
  \param iframe_max Actual rank last frame indice
  \param bolo_vect A vector containing the channel list (as a vector of string), for whole scan
  \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int write_data_flag_to_dirfile(struct param_common dir, struct samples samples_struct, long iframe_min, long iframe_max, std::vector<std::vector<std::string> > bolo_vect);

//! Reads RA and DEC tables in input scans and write them down to disk using getdata and dirfile format
/*!
  \param dir The param_common structure
  \param samples_struct The samples structure
  \param iframe_min Actual rank first frame indice
  \param iframe_max Actual rank last frame indice
  \param bolo_vect A vector containing the channel list (as a vector of string), for whole scan
  \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int write_RA_DEC_to_dirfile(struct param_common dir, struct samples samples_struct, long iframe_min, long iframe_max, std::vector<std::vector<std::string> > bolo_vect);

//! Reads data table in a dirfile pointed by "D" and stores it to "data" array
/*!
  \param D A pointer to an opened dirfile
  \param filename An input scan name (to find which data to open in the dirfile tree)
  \param field A channel name (to find which data to open in the dirfile tree)
  \param data An array containing the timeline from scan "filename" and channel "field"
  \param ns Scan "filename" number of samples
  \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_data_from_dirfile(DIRFILE* D, std::string filename, std::string field, double *&data, long ns);

//! Reads flag table in a dirfile pointed by "H" and stores it to "mask" array
/*!
  \param H A pointer to an opened dirfile
  \param filename An input scan name (to find which data to open in the dirfile tree)
  \param field A channel name (to find which data to open in the dirfile tree)
  \param mask An array containing the flags from scan "filename" and channel "field"
  \param ns Scan "filename" number of samples
  \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_flag_from_dirfile(DIRFILE* H, std::string filename, std::string field, int *&mask, long ns);

//! Reads RA table in a dirfile pointed by "D" and stores it to "ra" array
/*!
  \param D A pointer to an opened dirfile
  \param filename An input scan name (to find which data to open in the dirfile tree)
  \param field A channel name (to find which data to open in the dirfile tree)
  \param ra An array containing the RA coordinates from scan "filename" and channel "field"
  \param ns Scan "filename" number of samples
  \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_RA_from_dirfile(DIRFILE* D, std::string filename, std::string field, double *&ra, long ns);

//! Reads DEC table in a dirfile pointed by "D" and stores it to "dec" array
/*!
  \param D A pointer to an opened dirfile
  \param filename An input scan name (to find which data to open in the dirfile tree)
  \param field A channel name (to find which data to open in the dirfile tree)
  \param dec An array containing the DEC coordinates from scan "filename" and channel "field"
  \param ns Scan "filename" number of samples
  \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_DEC_from_dirfile(DIRFILE* D, std::string filename, std::string field, double *&dec, long ns);

//! Writes indpix to disk, in a binary file named indpix.bi
/*!
  \param ind_size Indpix size
  \param npix The number of filled pixels
  \param indpix Pixels indices in the map
  \param outdir The ouput directory pathname
  \param flagon Indicates whether a sample has been rejected or not
  \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int write_indpix(long long ind_size, long long npix, long long *indpix,  std::string outdir, int flagon);

//! Reads indpix from disk, in a binary file named indpix.bi
/*!
 \param ind_size Indpix size, read in the file
 \param npix The number of filled pixels
 \param indpix Pixels indices in the map, read in the file
 \param outdir The ouput directory pathname
 \param flagon Indicates whether a sample has been rejected or not, read in the file
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_indpix(long long &ind_size, long long &npix, long long *&indpix, std::string outdir, int &flagon);

//! Writes indpsrc to disk, in a binary file named "indpsrc.bi"
/*!
 \param map_size is indpsrc's size
 \param npixsrc The number of pixels in the mask
 \param indpsrc The masked pixels indices table
 \param outdir The ouput directory pathname
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int write_indpsrc(long long map_size,  long long  npixsrc, long long * indpsrc,   std::string outdir);

//! Reads indpsrc from disk, in a binary file named "indpsrc.bi"
/*!
 \param map_size is indpsrc's readed size, read from disk
 \param npixsrc The number of pixels in the mask, read from disk
 \param indpsrc The masked pixels indices table, to be read from disk
 \param outdir The ouput directory pathname
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int  read_indpsrc(long long &map_size, long long &npixsrc, long long *&indpsrc,   std::string outdir);

//! Writes PNd to disk, in a binary file named "filename"
/*!
 \param PNd is the Projected Noised data array
 \param npix The number of filled pixels
 \param filename The name of the binary file in which PNd will be written
 \param outdir The ouput directory pathname
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int write_PNd(double *PNd, long long npix, std::string outdir, std::string filename);

//! Reads PNd from disk, in a binary file named "filename"
/*!
 \param PNdtot is the Projected Noised data array
 \param npix The number of filled pixels
 \param filename The name of the binary file in which PNd will be written
 \param outdir The ouput directory pathname
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_PNd(double *&PNdtot, long long &npix, std::string outdir, std::string filename);

//! Writes samptopix to disk, in the dirfile pointed by "D"
/*!
 \param D A pointer to an opened dirfile
 \param ns The number of samples in "filename"
 \param samptopix Sample to pixel projection matrix, written to disk
 \param filename The considered scan file name (to find which data to open in the dirfile tree)
 \param boloname A channel name (to find which data to open in the dirfile tree)
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int write_samptopix(DIRFILE *D, long ns, long long *&samptopix, std::string filename, std::string boloname);

//! Reads samptopix from disk, in the dirfile pointed by "D"
/*!
 \param D A pointer to an opened dirfile
 \param ns The number of samples in "filename"
 \param samptopix Sample to pixel projection matrix, read from disk
 \param filename The considered scan file name (to find which data to open in the dirfile tree)
 \param boloname A channel name (to find which data to open in the dirfile tree)
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_samptopix(DIRFILE* D, long ns, long long *&samptopix, std::string filename, std::string boloname);

//! Writes a fourier transform to disk, in the dirfile pointed by "D"
/*!
 \param D A pointer to an opened dirfile
 \param ns The number of samples in "filename"
 \param fdata A fourier transform of the data, written in the disk
 \param prefixe A file prefix (fdata_ or fPs_) to allow 2 types of fourier transform data on disk at the same time
 \param filename The considered scan file name (to find which data to open in the dirfile tree)
 \param bolonames A channel list (to find which data to open in the dirfile tree)
 \param idet A channel index (to find which channel to consider in the list "bolonames")
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int write_fdata(DIRFILE *D, long ns, fftw_complex *fdata, std::string prefixe, long idet, std::string filename, std::vector<std::string> bolonames);

//! Reads a fourier transform from disk, in the dirfile pointed by "D"
/*!
 \param D A pointer to an opened dirfile
 \param ns The number of samples in "filename"
 \param fdata A fourier transform of the data, read from disk
 \param prefixe A file prefix (fdata_ or fPs_) to allow 2 types of fourier transform data on disk at the same time
 \param filename The considered scan file name (to find which data to open in the dirfile tree)
 \param bolonames A channel list (to find which data to open in the dirfile tree)
 \param idet A channel index (to find which channel to consider in the list "bolonames")
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_fdata(DIRFILE* D, long ns, fftw_complex *&fdata, std::string prefixe, long idet, std::string filename, std::vector<std::string> bolonames);

//! Reads a mixing matrix in a .txt file
/*!
 \param MixMatfile This file name contains the mixing matrix values
 \param ndet The number of detectors contained in the mixing matrix (determines mixing matrix size)
 \param ncomp The number of noise component to estimate by sanePS (determines mixing matrix size)
 \param mixmat The mixing matrix data (read from the file)
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_mixmat_txt(std::string MixMatfile, long ndet, long ncomp, double **&mixmat);



#endif
