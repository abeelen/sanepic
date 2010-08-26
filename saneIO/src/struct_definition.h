#ifndef STRUCT_DEFINITION_H_
#define STRUCT_DEFINITION_H_

#include <vector>
#include <string>

struct corner{
	double x;
	double y;
};

struct box {
	struct corner blc;
	struct corner trc;
};

struct common
/*! A structure that contains everything about directories, channel list and frame list */
{
	std::string dirfile;
	std::string output_dir;
	std::string tmp_dir;
	std::string noise_dir;

	std::string channel;
	std::string fits_filelist;
};

struct samples
/*! A structure that contains everything about frames, noise files and frame processing order */
{
	std::vector<std::string> fitsvect;
	std::vector<std::string> noisevect;
	std::vector<int> scans_index;

	std::vector<std::string> mixmat_file;
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

struct checksum
/*! A structure that contains sanePic input checksum values for crash recovery procedure */
{

	unsigned int chk_ini_file;
	unsigned int chk_wcs_file;
	unsigned int chk_pnd;
	unsigned int chk_indpix;
	unsigned int chk_indpsrc;

};

#endif /* STRUCT_DEFINITION_H_ */
