#ifndef MAP_MAKING_H_
#define MAP_MAKING_H_

//#include <cstdlib>
#include <string>
#include <vector>
#include "mpi_architecture_builder.h"


extern "C" {
#include <wcslib/wcs.h>
}


//! Computes map extrema by projecting the bolometers offsets back into the sky plane
/*!
 * outputs or updates (ra|dec)_(min|max) for SANEPIC format ONLY
 * Each rank treats his scans (indice_min and _max are given has inputs)
 \param samples_struct A samples structure
 \param iframe_min Actual rank first frame indice
 \param iframe_max Actual rank last frame indice
 \param ra_min ra minimum coordinate of the map : To be filled/replaced
 \param ra_max ra maximum coordinate of the map : To be filled/replaced
 \param dec_min dec minimum coordinate of the map : To be filled/replaced
 \param dec_max dec maximum coordinate of the map : To be filled/replaced
 \param bolo_vect A vector containing the channel list (as a vector of string), for whole scan
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int computeMapMinima(struct samples samples_struct, std::string dirfile,
		long iframe_min, long iframe_max,
		double &ra_min,double &ra_max,double &dec_min,double &dec_max, std::vector<std::vector<std::string> > bolo_vect);

//! Computes map extrema by projecting the bolometers offsets back into the sky plane
/*!
 * outputs or updates (ra|dec)_(min|max) for HIPE format ONLY
 * Each rank treats his scans (indice_min and _max are given has inputs)
 \param tmp_dir A string containing the temporary files pathname
 \param samples_struct A samples structure
 \param iframe_min Actual rank first frame indice
 \param iframe_max Actual rank last frame indice
 \param ra_min ra minimum coordinate of the map : To be filled/replaced
 \param ra_max ra maximum coordinate of the map : To be filled/replaced
 \param dec_min dec minimum coordinate of the map : To be filled/replaced
 \param dec_max dec maximum coordinate of the map : To be filled/replaced
 \param bolo_vect A vector containing the channel list (as a vector of string), for whole scan
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int computeMapMinima_HIPE(std::string tmp_dir, struct samples samples_struct,
		long iframe_min, long iframe_max,
		double &ra_min,double &ra_max,double &dec_min,double &dec_max, std::vector<std::vector<std::string> > bolo_vect);


//! Computes an array's minimum and maximum values without taking flagged samples into account
/*!
 * Fills min_array and max_array
 \param array An array containing RA or DEC coordinates
 \param flag An array containing the flags corresponding to "array" (same channel)
 \param size Array's size (flag's size is the same)
 \param min_array Array's Minimum value is stored in this variable
 \param max_array Array's Maximum value is stored in this variable
 \return An integer >0 if the whole flag values are >0 (array is fully flagged), and 0 otherwise
 */
int minmax_flag(double  *& array, int *& flag, long size, double & min_array, double &  max_array);


//! Computes Map Header using projection and map informations
/*!
 * Fills wcsprm struct
 \param pixdeg The pixel size in the map
 \param ctype A char* specifying axis labels, or coordinates type (default is "EQ" : equatorial system RA/DEC)
 \param prjcode A char* specifying projection type (default is "TAN" : tangent plane)
 \param coordscorner An array containing ra/dec min/max of the final map
 \param wcs A pointer to a wcsprm struct that will be filled in this routine
 \param NAXIS1 Number of horizontal pixels
 \param NAXIS2 Number of vertical pixels
 */
void computeMapHeader(double pixdeg, char *ctype, char* prjcode, double * coordscorner,
		struct wcsprm *& wcs, long &NAXIS1, long &NAXIS2);


//! Computes a Naïve map using preprocessing and a projection of the raw data
/*!
 * Computes Naïve map, coverage map and a history table (also copy binary mask if any)
 \param samples_struct A samples structure
 \param PNd Projected Noised data
 \param outdir A string containing the output pathname
 \param det A vector containing the channel list (as a vector of string), for actual scan (referenced by iframe)
 \param ndet The number of channels contained in det
 \param orderpoly Polynomia order to be removed from the timeline
 \param napod Number of samples to apodize at data begin and end
 \param f_lppix High-pass Filter cut-off frequency (converted in samples)
 \param ns Actual frame's (files[iframe]) number of samples
 \param indpix The pixels indices table
 \param iframe Actual frame loop indice
 \param hits This array corresponds to map's coverage. Stored in naïve map has a fits image
 \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int do_PtNd_Naiv(struct samples samples_struct, double *PNd, std::string outdir, std::vector<std::string> det, long ndet, int orderpoly, int napod, double f_lppix, long ns,
		long long *indpix, long iframe, long *hits);

#endif /* MAP_MAKING_H_ */
