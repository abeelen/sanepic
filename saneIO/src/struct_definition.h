#ifndef STRUCT_DEFINITION_H_
#define STRUCT_DEFINITION_H_

#include <vector>
#include <string>

struct common
/*! A structure that contains everything about directories, channel list and frame list */
{
	std::string dirfile;
	std::string output_dir;
	std::string tmp_dir;
	std::string noise_dir;
	std::string input_dir;

	std::string channel;
	std::string fits_filelist;
	std::string bolo_global_filename;
	std::string suffix;
};

struct samples
/*! A structure that contains everything about frames, noise files and frame processing order */
{
	std::vector<std::string> fitsvect;
	std::vector<std::string> noisevect;
	std::vector<int> scans_index;

	std::string cov_matrix_file;

	bool framegiven;
	std::string *fits_table;
	std::string *noise_table;
	int *index_table;
	long *nsamples;
	long ntotscan;
	std::string filename; // What is this ? : name of the fits_filelist.txt file read in ini file !
};


struct detectors
/*! A structure that contains the name of the detectors + number of det */
{
	long ndet;
	std::vector<std::string> boloname;
};

struct param_positions
/*! A structure that contains user options about map projection and properties */
{
	std::string maskfile;
	double pixdeg;
	bool flgdupl;
	bool projgaps;
	int fileFormat;
};

struct param_process
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

struct PS
/*! A structure that contains user options about sanePS procedure */
{
	double fcutPS;
	std::string ell_suffix;
	std::string mix_suffix;
	std::string ell_global_file;
	std::string mix_global_file;
	std::vector<std::string> ell_names;
	std::vector<std::string> mix_names;
	std::string signame;
	long ncomp;
};

struct sanePic
/*! A structure that contains user options about sanePic procedure */
{
	int iterw;
	int itermax;
	int save_data;
	int restore;
	double thresholds; // determine thresholds
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
