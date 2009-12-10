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

struct input_commons {
	long napod;
	bool NOFILLGAP;
	bool flgdupl;
	double pixdeg;
};

struct user_options {
	double fsamp;
	bool NORMLIN;
	bool remove_polynomia;
	int poly_order;
	bool CORRon;
	double f_lp;
	bool projgaps;
};


#endif /* STRUCT_DEFINITION_H_ */
