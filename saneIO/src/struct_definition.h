#ifndef STRUCT_DEFINITION_H_
#define STRUCT_DEFINITION_H_


struct corner{
	double x;
	double y;
};

struct box {
	struct corner blc;
	struct corner trc;
};

struct common {
	std::string dirfile;
	std::string output_dir;
	std::string tmp_dir;
	std::string noise_dir;

	std::string channel;
	std::string fits_filelist;
};

struct samples {
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
	std::string filename; //TODO: What is this ?
};

struct detectors {
	long ndet;
	std::vector<std::string> boloname;
};

struct param_positions {
	std::string maskfile;
	double pixdeg;
	bool flgdupl;
	bool projgaps;
	int fileFormat;
};

struct param_process {
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


#endif /* STRUCT_DEFINITION_H_ */
