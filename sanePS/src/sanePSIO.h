/*
 * sanePSIO.h
 *
 *  Created on: Dec 27, 2012
 *      Author: abeelen
 */

#ifndef SANEPSIO_H_
#define SANEPSIO_H_

#include <vector>
#include <string>

#include "StructDefinition.h"

//! restore a sanePS session using temporary session state values (completed_step, commonm2, N, P, Rellth, Rellexp, SPref) depending on which step has to be restore
/*!
\param tmp_dir The temporary files directory
\param filename The Scan that was estimated before the session crashed
\param completed_step The number of sanePS steps that were completed before the crash
\param commonm2 sanePS internal Matrix data
\param N sanePS internal Matrix data
\param P sanePS internal Matrix data
\param Rellth sanePS internal Matrix data (Theorical noise)
\param Rellexp sanePS internal Matrix data (Experimental noise)
\param SPref sanePS internal array data
\param ndet number of channels for the considered scan "filename"
\param ncomp number of noise component to estimate
\param nbins number of spectra bins
\param ns Number of samples in the scan named "filename"
\return An integer >0 if there were a problem, or 0 if everything went OK
*/
int restore_session(std::string tmp_dir, std::string filename, int &completed_step, double **commonm2, double **N,
		double **P, double **Rellth, double **Rellexp, double *SPref, long ndet, int ncomp, long nbins, long ns);

//! Save a sanePS session using temporary session state values (completed_step, commonm2, N, P, Rellth, Rellexp, SPref) depending on which step has to be restore
/*!
\param tmp_dir The temporary files directory
\param filename The Scan that was estimated before the session crashed
\param completed_step The number of sanePS steps that were completed before the crash
\param commonm2 sanePS internal Matrix data
\param N sanePS internal Matrix data
\param P sanePS internal Matrix data
\param Rellth sanePS internal Matrix data (Theorical noise)
\param Rellexp sanePS internal Matrix data (Experimental noise)
\param SPref sanePS internal array data
\param ndet number of channels for the considered scan "filename"
\param ncomp number of noise component to estimate
\param nbins number of spectra bins
\param ns Number of samples in the scan named "filename"
\return An integer >0 if there were a problem, or 0 if everything went OK
*/
int save_session(std::string tmp_dir, std::string filename, int completed_step, double **commonm2, double **N,
		double **P, double **Rellth, double **Rellexp, double *SPref, long ndet, int ncomp, long nbins, long ns);


int write_to_disk(std::string outdirSpN, std::string fits_filename, struct param_sanePS structPS, std::vector<std::string> det, long nbins, double *ell, double **mixmat,
		double **Rellth, double **Rellexp, double **N, double *SPref, double **P);

int write_MixMatrix(std::string fname, std::vector<std::string> bolos, long ncomp, double **mixmat);

int read_MixMatrix(string fname, std::vector<string> &det_vect, long &ncomp, double **& mixmat);

uint16_t assignMixMat(std::string fname, std::vector<std::string> det, long ncomp, double **& mixmat);

//! Used in sanePS to write DEBUG channels psd
/*!
 \param fname output fits filename
 \param nx Number of lines (fits image)
 \param ny Number of columns (fits image)
 \param dtype 'd' for double, 'l' for long
 \param psd1d The array to write in the fits image
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */

int write_psd_tofits(std::string fname, long nx, long ny, char dtype, void * psd1d, std::string ext_name, bool extend);


//! Reads a mixing matrix in a .txt file
/*!
 \param MixMatfile This file name contains the mixing matrix values
 \param ndet The number of detectors contained in the mixing matrix (determines mixing matrix size)
 \param ncomp The number of noise component to estimate by sanePS (determines mixing matrix size)
 \param mixmat The mixing matrix data (read from the file)
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
//int readMixmatTxt(std::string MixMatfile, long ndet, long ncomp, double **&mixmat);


#endif /* SANEPSIO_H_ */
