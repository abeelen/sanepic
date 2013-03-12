#ifndef STRUCT_DEFINITION_H_
#define STRUCT_DEFINITION_H_

#include <vector>
#include <string>

#include <fftw3.h>


extern "C"{
#include "getdata.h"
}

struct param_common
/*! A structure that contains everything about directories, channel list and frame list */
{
	std::string data_dir; /*! data directory (fits files) */
	std::string output_dir; /*! output directory */
	std::string tmp_dir; /*! temporary directory */
	std::string input_dir; /*! input directory (config files) */

	std::string fits_filelist; /*! A list (txt file) in which are listed fits files (and processors orders if MPI : optional) */

	std::string bolo; /*! A list (txt file) in which are listed used channels for the whole scans (optional) */
	std::string bolo_suffix;  /*! A string suffix that indicates what to add to scan names to find scans channel list (one list per scan, optional) */

	std::string parallel_scheme; /*! A string defining how to parallelize things (2 possibles values: bolo, frame)*/
};

struct param_sanePos
/*! A structure that contains user options about map projection and properties */
{
	std::string maskfile; /*! A fits file in which a binary mask of the map is given */

	double pixdeg; /*! The pixel size */
	double lon; /*! Ra nominal (optional) */
	double lat;  /*! Dec nominal (optional) */

	std::string radesys; /*! Coordinate reference frame for the RA and DEC */
	double equinox; /*! Equinox of celestial coordinate system */
	double restwav; /*! rest wavelength in vacuo [m]" */

	std::string axistype;     /*! Axis type (EQ/Gal, default is EQ) */
	std::string projcode; /*! Projection code (Default is TAN, optional) */

	bool eq2gal; /* Projection from J2000.0 RA/DEC to Galactic coordinates (default False) */
	bool gal2eq; /* Projection from Galactic Coordinates to J2000.0 RA/DEC (default False) */

	bool flgdupl; /*! True if flagged data are put in a separated map (default is False) */
	bool projgaps; /*! Gaps are projected in the maps ? */
	int fileFormat; /*! indicates fits files format : HIPE = 1, SANEPIC = 0 */
};

struct param_saneProc
/*! A structure that contains user options about preprocessing properties */
{
	bool remove_linear ;   /*! True if a simple baseline removed from the data */
	bool fill_gap;         /*! True if data gaps are filled */
	bool CORRon;           /*! True if correlation between detectors are included in the analysis */
	bool remove_polynomia; /*! True if a polynomia is removed from the data : poly_order > 0 */
	bool highpass_filter;  /*! True if an high pass filter is perfomed on the data */

	long napod;            /*! number of samples to apodize at data begin and end */
	int poly_order;        /*! Polynomia order to be removed from the timeline */

	double fsamp;          /*! sampling frequency (Hz) (same for all scans) */
	double fhp;            /*! high pass frequency filter for the time stream (same for all scans)*/
	double fcut;           /*! noise filter cut-off frequency (same for all scans) */

	std::string fsamp_file; /*! sampling frequencies (per scan) (Hz) */
	std::string fhp_file;   /*! high pass frequencies filter for the timestream (per scan) */
	std::string fcut_file;  /*! noise filter cut-off frequencies (per scan) */


	bool wisdom ;          /*! True will force sanePre to precompute some fftw to acquire wisdom */


};

struct param_saneInv
/*! A structure that contains informations about covariance matrices filenames */
{
	std::string noise_dir; /*! covariance matrices directory */

	std::string cov_matrix;  /*! covariance matrix filename (same for all scans) */
	std::string cov_matrix_suffix; /*! add this suffix to fits filenames to obtain the corresponding cov matrix */
};

struct param_sanePS
/*! A structure that contains user options about sanePS procedure */
{
	std::string ell_suffix; /*! suffix to add to fits filenames to obtain the bins for the noise spectrum */
	std::string mix_suffix; /*! suffix to add to fits filenames to obtain mixing matrix for the given scan */
	std::string ell; /*! file containing the bins for the noise spectrum (same for all scans) */
	std::string mix; /*! file containing the mixing matrix (same for all scans) */

	//TODO: Ugly turnaround until sanePS is released;
	std::string cov_matrix;  /*! covariance matrix filename (same for all scans) */
	std::string cov_matrix_suffix;  /*! add this suffix to fits filenames to obtain the corresponding cov matrix */


	std::string signame; /*! fits file containing the map that should be substracted to the data for a better Noise estimation */
	int niter; /*! number of iteration for the expectation minimization step */
	int ncomp; /*! number of component(s) to estimate */
	bool restore; /*! If true, restore a previous session for sanePS */
	bool save_data; /*! If true, save the actual session for sanePS after each step */
};

struct param_sanePic
/*! A structure that contains user options about sanePic procedure */
{
	int iterw; /*! Write temporary map files on disk every iterW number of loop */
	int itermax; /*! Maximum iteration number for conjugate gradient */
	double tolerance; /*! Tolerance to reach to stop iterations */
	double subtolerance; /*! Tolerance for the first iteration */
	int save_data; /*! If true, save the actual session for sanePic after each iterw loop */
	int restore; /*! If true, restore a previous session for sanePic */
	std::string map_prefix;  /*! add this prefix to the maps generated by sanePic */
	//	double thresholds; // determine thresholds
};

struct samples
/*! A structure that contains everything about frames, noise files and frame processing order */
{
	bool framegiven;                    /*!  True if the processor/scan index is given in fitsfilelist */

	long iframe_min, iframe_max;        /*! indexes  to process in the following vector */

	std::vector<DIRFILE *> dirfile_pointers; /*! pointers to the dirfiles that contains temporary binaries */

	std::vector<long> scans_index;       /*!  a vector containing the index of the scans */

	std::vector<std::string> fitsvect;  /*! a vector containing input fits filenames */
	std::vector<std::string> noisevect; /*! a vector containing input covariance matrices filenames */
	std::vector<std::string> basevect;  /*! a vector containing input fits basenames : "." are replaced by _ and ".fits" is removed */
	std::vector<std::string> bolovect;  /*! a vector containing input bolometer lists filenames */
	std::vector<std::vector<std::string> > bolo_list;

	std::vector<double> fcut;           /*! a vector containing Noise filter cut-off frequency for each scan */
	std::vector<double> fsamp;          /*! a vector containing the sampling frequencies for each scan */
	std::vector<double> fhp;            /*! a vector containing data high pass filter frequency for each scan */

	std::vector<long> nsamples;         /*! a vector containing the number of samples for each input fits filenames */
	long ntotscan;                      /*! the total number of scans, should be the size of all vector in this struct */

	// for sanePS only
	std::vector<std::string> ell_names; /*! a vector for the bins (for the noise spectrum) filenames */
	std::vector<std::string> mix_names; /*! a vector containing the mixing matrices filenames */
	std::vector<long> nbins;            /*!  a vector containing the number of ell bins for each scan */
	std::vector<long> ndet;             /*! a vector containing the number of detector for each scan */

	uint parallel_scheme;

};

struct checksum
/*! A structure that contains sanePic input checksum values for crash recovery procedure */
{

	unsigned int chk_ini_file; /*! the whole structure is used to compute this checksum */
	unsigned int chk_wcs_file; /*! mapheader.keyrec checksum */
	unsigned int chk_pnd; /*! PNd.bi checksum */
	unsigned int chk_indpix; /*! indpix.bi checksum */
	unsigned int chk_indpsrc; /*! indpsrc.bi checksum */

};

struct checkHDU
/*! A structure that determines which tables saneCheck has to check */
{
	bool checkREFERENCEPOSITION; /*! True if saneCheck has to check reference position table */
	bool checkOFFSETS; /*! True if saneCheck has to check offsets table */
	bool checkLON; /*! True if saneCheck has to check RA table */
	bool checkLAT; /*! True if saneCheck has to check DEC table */
};

struct param_saneCheck
/*! A structure that contains saneCheck specific ini file informations */
{
	struct checkHDU Check_it; /*! A structure that determines which tables saneCheck has to check */
	std::vector<double>  bolo_gain_check;  /*! A vector that contains bolometers gains */
	bool checkNAN; /*! True if saneCheck has to check NANs in each table */
	bool checktime; /*! True if saneCheck has to check time gaps in time table */
	bool checkGain;  /*! True if saneCheck has to recompute bolometers gains */
	bool checkflag; /*! True if saneCheck has to check time for fully flagged detectors */
};

#endif /* STRUCT_DEFINITION_H_ */
