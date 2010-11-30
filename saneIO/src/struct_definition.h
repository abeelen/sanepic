#ifndef STRUCT_DEFINITION_H_
#define STRUCT_DEFINITION_H_

#include <vector>
#include <string>

struct param_common
/*! A structure that contains everything about directories, channel list and frame list */
{
	std::string dirfile;
	std::string output_dir;
	std::string tmp_dir;
	std::string noise_dir;
	std::string input_dir;

	std::string fits_filelist;
	std::string bolo_global_filename;
	std::string bolo_suffix;
};

struct param_sanePos
/*! A structure that contains user options about map projection and properties */
{
	std::string maskfile;

	double pixdeg;
	double ra_nom;
	double dec_nom;
	std::string projtype;

	bool flgdupl;
	bool projgaps;
	int fileFormat;
};

struct param_sanePre
/*! A structure that contains user options about preprocessing properties */
{
	bool NORMLIN;
	bool NOFILLGAP;
	bool CORRon;
	bool remove_polynomia;

	long napod;
	int poly_order;
	double fsamp;
	double f_lp;

	std::string fcut_file;
};

struct param_saneInv
/*! A structure that contains informations about covariance matrices filenames */
{
	std::string cov_matrix_file;
	std::string cov_matrix_suffix;
};

struct param_sanePS
/*! A structure that contains user options about sanePS procedure */
{
	std::string ell_suffix;
	std::string mix_suffix;
	std::string ell_global_file;
	std::string mix_global_file;

	//TODO: Ugly turnaround until sanePS is released;
	std::string cov_matrix_file;
	std::string cov_matrix_suffix;


	std::string signame;
	int ncomp;
};

struct param_sanePic
/*! A structure that contains user options about sanePic procedure */
{
	int iterw;
//	int itermax; // TODO : add itermax + thresholds in sanepic_ini
	int save_data;
	int restore;
//	double thresholds; // determine thresholds
};

struct samples
/*! A structure that contains everything about frames, noise files and frame processing order */
{
	std::vector<std::string> fitsvect;
	std::vector<std::string> noisevect;
	std::vector<std::string> bolovect;


	std::vector<double> fcut;

	std::vector<std::string> ell_names;
	std::vector<std::string> mix_names;

	std::vector<int> scans_index;

	bool framegiven;

	//TODO: Switch nsample to vector too...
	std::vector<long> nsamples;
	long ntotscan;
};

struct checksum
/*! A structure that contains sanePic input checksum values for crash recovery procedure */
{

	unsigned int chk_ini_file;
	unsigned int chk_wcs_file;
	unsigned int chk_pnd;
	unsigned int chk_indpix;
	unsigned int chk_indpsrc;

};

struct checkHDU {
	bool checkREFERENCEPOSITION;
	bool checkOFFSETS;
	bool checkRA;
	bool checkDEC;
};

struct saneCheck
/*! A structure that contains saneCheck specific ini file informations */
{
	struct checkHDU Check_it;
	std::vector<double>  bolo_gain_check;
	bool checkNAN;
	bool checktime;
	bool checkGain;
	bool checkflag;


};

#endif /* STRUCT_DEFINITION_H_ */
