#ifndef CORR_PREPROCESS_H_
#define CORR_PREPROCESS_H_

#include <string>
#include <fstream>
#include <vector>
#include "struct_definition.h"
#include <fftw3.h>

extern "C" {
#include "nrutil.h"
}

//! Computes, and stores to disk, Fourier transform of the data after pre-processing
/*!
 * S is deprojected and removed from the data
 \param S Map signal
 \param tmp_dir A string containing the temporary files pathname
 \param dirfile fits files data directory
 \param samples_struct A samples structure
 \param pos_param The param_sanePos structure
 \param proc_param The param_sanePre structure
 \param NAXIS1 Number of horizontal pixels (determined by sanePos)
 \param NAXIS2 Number of vertical pixels (determined by sanePos)
 \param addnpix Number of pix to add to compute the final maps in case of duplication and/or user binary mask
 \param indpsrc The masked pixels indices table, to be read from disk
 \param npixsrc number of pixels included in Crossing Constraint Removal (binary mask)
 \param det A vector containing the channel list (as a vector of string), for whole scan
 \param ndet The number of detector in det
 \param indpix The pixels indices table
 \param npix Number of filled pixels (in the map)
 \param f_lppix High-pass Filter cut-off frequency (converted in samples)
 \param ns Number of samples for the considered scan : samples_struct.fitsvect[iframe]
 \param iframe Current loop frame indice
 \param para_bolo_indice In case para_bolo is used : processor rank is given to select which channels it has to computes
 \param para_bolo_size In case para_bolo is used : processor size is given to computes channels indexes that has to be computed
 \param fname log file that is used in debug mode to track processing time between detectors
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int write_ftrProcesdata(double *S, struct param_sanePre proc_param, struct samples samples_struct, struct param_sanePos pos_param,
		std::string tmp_dir, std::vector<std::string> det, long ndet, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2,
		long long npix,	long long npixsrc, long long addnpix, double f_lppix, long ns, long iframe, int para_bolo_indice, int para_bolo_size, std::string fname);

//! Computes Mp and hits arrays using write_ftrProcesdata or write_tfAS outputs
/*!
 \param samples_struct A samples structure
 \param PNd is the Projected Noised data array
 \param prefixe A string prefixe that is used to determine which fourier transform binary file to open (fdata or fPs)
 \param NAXIS1 Number of horizontal pixels (determined by sanePos)
 \param NAXIS2 Number of vertical pixels (determined by sanePos)
 \param det A vector containing the channel list (as a vector of string), for whole scan
 \param ndet The number of detector in det
 \param indpix The pixels indices table
 \param npix Number of filled pixels (in the map)
 \param f_lppix High-pass Filter cut-off frequency (converted in samples)
 \param fsamp Instrument Sampling frequency
 \param ns Number of samples for the considered scan : samples_struct.fitsvect[iframe]
 \param iframe Current loop frame indice
 \param para_bolo_indice In case para_bolo is used : processor rank is given to select which channels it has to computes
 \param para_bolo_size In case para_bolo is used : processor size is given to computes channels indexes that has to be computed
 \param Mp CG preconditionner, outputed by do_PtNd
 \param hits Map coverage, outputed by do_PtNd
 \param fname log file that is used in debug mode to track processing time between detectors
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int do_PtNd(struct samples samples_struct, double *PNd, std::string prefixe,
		std::vector<std::string> det, long ndet, double f_lppix, double fsamp, long ns, int para_bolo_indice, int para_bolo_size,
		long long *indpix, long NAXIS1, long NAXIS2, long long npix, long iframe,
		double *Mp, long *hits,std::string fname);

//! Computes, and stores to disk, Fourier transform of the signal S after deprojection
/*!
 \param S Map signal
 \param samples_struct A samples structure
 \param NAXIS1 Number of horizontal pixels (determined by sanePos)
 \param NAXIS2 Number of vertical pixels (determined by sanePos)
 \param det A vector containing the channel list (as a vector of string), for whole scan
 \param ndet The number of detector in det
 \param indpix The pixels indices table
 \param npix Number of filled pixels (in the map)
 \param flgdupl True if flagged data are put in a separated map
 \param ns Number of samples for the considered scan : samples_struct.fitsvect[iframe]
 \param filename Current loop fits filename : samples_struct.fitsvect[iframe]
 \param para_bolo_indice In case para_bolo is used : processor rank is given to select which channels it has to computes
 \param para_bolo_size In case para_bolo is used : processor size is given to computes channels indexes that has to be computed
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int write_tfAS(struct samples samples_struct, double *S, std::vector<std::string> det, long ndet, long long *indpix, long NAXIS1, long NAXIS2, long long npix,
		bool flgdupl, long ns, std::string filename, int para_bolo_indice, int para_bolo_size);

#endif /* CORR_PREPROCESS_H_ */
