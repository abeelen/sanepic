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


struct directories {
	std::string dirfile;
	std::string outdir;
	std::string tmp_dir;
};

struct samples {
	std::vector<std::string> fitsvect;
	std::vector<std::string> noisevect;
	std::vector<long> scans_index;
	bool framegiven;
	std::string *fits_table;
	std::string *noise_table;
	long *index_table;
	long *nsamples;
	long ntotscan;
	std::string filename;
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
};

struct param_process {
	long napod;
	bool NOFILLGAP;
	double fsamp;
	bool NORMLIN;
	bool remove_polynomia;
	int poly_order;
	bool CORRon;
	double f_lp;
};


#endif /* STRUCT_DEFINITION_H_ */
