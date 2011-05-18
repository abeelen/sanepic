#ifndef STRUCT_DEFINITION_H_
#define STRUCT_DEFINITION_H_

#include <vector>
#include <string>


extern "C"{
#include "getdata.h"
}

#define sanepic_version "0.5"

#ifdef USE_MPI

struct ini_var_strings
/*! A structure that contains strings sizes about project structures filled with ini file variables */
{
	// common
	int dirfile; /*! common.dirfile.size() */
	int output_dir; /*! common.output_dir.size() */
	int input_dir; /*! common.input_dir.size() */
	int tmp_dir; /*! common.tmp_dir.size() */

	int fits_filelist; /*! common.fits_filelist.size() */
	int bolo_global_filename; /*! common.bolo_global_filename.size() */
	int bolo_suffix; /*! common.bolo_suffix.size() */

	// sanePos
	int maskfile; /*! param_position.maskfile.size() */
	int projtype; /*! param_position.projtype.size() */

	// sanePre
	int fcut_file;  /*! proc_param.fcut_file.size() */

	// saneInv
	int noise_dir; /*! saneInv_param.noise_dir.size() */
	int cov_matrix_file; /*! saneInv_param.cov_matrix_file.size() */
	int cov_matrix_suffix; /*! saneInv_param.cov_matrix_suffix.size() */

	// sanePS
	int ell_suffix;  /*! sanePS_param.ell_suffix.size() */
	int ell_global_file;  /*! sanePS_param.ell_global_file.size() */
	int signame; /*! sanePS_param.signame.size() */
	int mix_global_file; /*! sanePS_param.mix_global_file.size() */
	int mix_suffix; /*! sanePS_param.mix_suffix.size() */

	// samples
	int ntotscan; /*! samples_struct.ntotscan */
	int *fitsvect; /*! A table containing fitsvect strings sizes */
	int *noisevect;  /*! A table containing noisevect strings sizes */
	int *bolovect;  /*! A table containing bolovect strings sizes */


	// all projects
	int sizemax; /*! considering all those sizes, sizemax is the maximum */

};

#endif


struct param_common
/*! A structure that contains everything about directories, channel list and frame list */
{
	std::string dirfile; /*! data directory (fits files) */
	std::string output_dir; /*! output directory */
	std::string tmp_dir; /*! temporary directory */
	std::string input_dir; /*! input directory (config files) */

	std::string fits_filelist; /*! A list (txt file) in which are listed fits files (and processors orders if MPI : optional) */
	std::string bolo_global_filename; /*! A list (txt file) in which are listed used channels for the whole scans (optional) */
	std::string bolo_suffix;  /*! A string suffix that indicates what to add to scan names to find scans channel list (one list per scan, optional) */

};

struct param_sanePos
/*! A structure that contains user options about map projection and properties */
{
	std::string maskfile; /*! A fits file in which a binary mask of the map is given */

	double pixdeg; /*! The pixel size */
	double ra_nom; /*! Ra nominal (optional) */
	double dec_nom;  /*! Dec nominal (optional) */
	std::string projtype; /*! Projection type (Default is TAN, optional) */

	bool flgdupl; /*! True if flagged data are put in a separated map (default is False) */
	bool projgaps; /*! Gaps are projected in the maps ? */
	int fileFormat; /*! indicates fits files format : HIPE = 1, SANEPIC = 0 */
};

struct param_sanePre
/*! A structure that contains user options about preprocessing properties */
{
	bool NORMLIN; /*! True if a simple baseline removed from the data */
	bool NOFILLGAP; /*! True if data gaps are filled */
	bool CORRon; /*! True if correlation between detectors are included in the analysis */
	bool remove_polynomia; /*! True if a polynomia is removed from the data : poly_order > 0 */

	long napod; /*! number of samples to apodize at data begin and end */
	int poly_order; /*! Polynomia order to be removed from the timeline */
	double fsamp; /*! Detectors sampling frequency (Hz) */
	double f_lp; /*! frequency of the high pass filter applied to the data */

	std::string fcut_file; /*! Noise filter cut-off frequency */
};

struct param_saneInv
/*! A structure that contains informations about covariance matrices filenames */
{
	std::string noise_dir; /*! covariance matrices directory */

	std::string cov_matrix_file;  /*! covariance matrix filename (same for all scans) */
	std::string cov_matrix_suffix; /*! add this suffix to fits filenames to obtain the corresponding cov matrix */
};

struct param_sanePS
/*! A structure that contains user options about sanePS procedure */
{
	std::string ell_suffix; /*! suffix to add to fits filenames to obtain the bins for the noise spectrum */
	std::string mix_suffix; /*! suffix to add to fits filenames to obtain mixing matrix for the given scan */
	std::string ell_global_file; /*! file containing the bins for the noise spectrum (same for all scans) */
	std::string mix_global_file; /*! file containing the mixing matrix (same for all scans) */

	//TODO: Ugly turnaround until sanePS is released;
	std::string cov_matrix_file;  /*! covariance matrix filename (same for all scans) */
	std::string cov_matrix_suffix;  /*! add this suffix to fits filenames to obtain the corresponding cov matrix */


	std::string signame; /*! fits file containing the map that should be substracted to the data for a better Noise estimation */
	int ncomp; /*! number of component(s) to estimate */
	bool restore; /*! If true, restore a previous session for sanePS */
	bool save_data; /*! If true, save the actual session for sanePS after each step */
};

struct param_sanePic
/*! A structure that contains user options about sanePic procedure */
{
	int iterw; /*! Write temporary map files on disk every iterW number of loop */
	int itermax; /*! Maximum iteration number for conjugate gradient */
	int save_data; /*! If true, save the actual session for sanePic after each iterw loop */
	int restore; /*! If true, restore a previous session for sanePic */
	std::string map_prefix;  /*! add this prefix to the maps generated by sanePic */
	//	double thresholds; // determine thresholds
};

struct samples
/*! A structure that contains everything about frames, noise files and frame processing order */
{
	std::vector<std::string> fitsvect;  /*! a vector containing input fits filenames */
	std::vector<std::string> noisevect; /*! a vector containing input covariance matrices filenames */
	std::vector<std::string> bolovect; /*! a vector containing input bolometer lists filenames */
	std::vector<std::string> basevect; /*! a vector containing input fits basenames : "." are replaced by _ and ".fits" is removed */


	DIRFILE *dirfile_pointer; /*! a pointer to the dirfile that contains temporary binaries */

	std::vector<double> fcut; /*! a vector containing Noise filter cut-off frequency for each scan */

	std::vector<std::string> ell_names; /*! a vector for the bins (for the noise spectrum) filenames */
	std::vector<std::string> mix_names; /*! a vector containing the mixing matrices filenames */

	//noise binary sizes
	std::vector<long> nbins; /*!  a vector containing the number of bins for each scan */
	std::vector<long> ndet; /*! a vector containing the number of detector for each scan */

	std::vector<int> scans_index; /*!  a vector containing the index of the scans */

	bool framegiven; /*!  True if the processor/scan index is given in fitsfilelist */

	std::vector<long> nsamples; /*! a vector containing the number of samples for each input fits filenames */
	long ntotscan; /*! the total number of scans */
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
	bool checkRA; /*! True if saneCheck has to check RA table */
	bool checkDEC; /*! True if saneCheck has to check DEC table */
};

struct saneCheck
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
