
#ifndef COVMATRIX_IO_H_
#define COVMATRIX_IO_H_

#include <string>
#include <vector>

extern "C" {
#include "getdata.h"
}

//!  Writes the Inverse Covariance Matrices (and bins : ell) to disk, using a binary format
/*!
 \param D A pointer to an opened dirfile
 \param bolos A channel list (to write the data in the dirfile tree)
 \param nbins Number of spectra bins
 \param ell The spectra bins, written to disk
 \param Rellth Inversed Theorical Covariance matrix, written to disk
 \param suffix A file suffix to add to fits file name (_InvNoisePS : arbitrary chosen by us)
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int write_InvNoisePowerSpectra(DIRFILE* D, std::vector<std::string> bolos, long nbins, double * ell, double **Rellth, std::string suffix);

//!  Reads the Inverse Covariance Matrices (and bins : ell) from disk, using a binary format
/*!
 \param D A pointer to an opened dirfile
 \param boloName A channel name (to find which data to open in the dirfile tree)
 \param nbins Number of spectra bins
 \param ndet Number of channels to be used for the considered scan
 \param ell The spectra bins, read from disk
 \param SpN_all Inversed Theorical Covariance matrix, read from disk
 \param suffix A file suffix to add to fits file name (_InvNoisePS : arbitrary chosen by us)
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_InvNoisePowerSpectra(DIRFILE* D, std::string boloName, std::string suffix, long nbins, long ndet, double ** ell, double *** SpN_all);

//! This function writes the NoiseNoise Matrices in a fits file (also Writes ell)
/*!
 \param fname output fits filename
 \param bolos A channel list
 \param nbins Number of spectra bins
 \param ell The spectra bins, written to fits file
 \param Rellth Theorical Covariance matrix, written to fits file
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int write_CovMatrix(std::string fname, std::vector<std::string> bolos, long nbins, double *ell, double **Rellth);

//! This function reads the NoiseNoise Matrices in a fits file (also reads the mixing matrices)
/*!
 \param fname Covariance Matrix fits filename
 \param bolos A channel list (filled by read_CovMatrix)
 \param nbins Number of spectra bins (filled by read_CovMatrix)
 \param ell The spectra bins, read in fits file (filled by read_CovMatrix)
 \param Rellth Theorical Covariance matrix, read in fits file (filled by read_CovMatrix)
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int read_CovMatrix(std::string fname, std::vector<std::string> &bolos, long &nbins, double *&ell, double **&Rellth);

//! Used in sanePS to write DEBUG channels psd
/*!
 \param fname output fits filename
 \param nx Number of lines (fits image)
 \param ny Number of columns (fits image)
 \param dtype 'd' for double, 'l' for long
 \param psd1d The array to write in the fits image
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int write_psd_tofits(std::string fname, long nx, long ny, char dtype, void * psd1d);


#endif /* COVMATRIX_IO_H_ */
