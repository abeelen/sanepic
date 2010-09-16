

#ifndef TOOLS_H_
#define TOOLS_H_

struct checkHDU {
	bool checkREFERENCEPOSITION;
	bool checkOFFSETS;
	bool checkRA;
	bool checkDEC;
};

/*! this function determines which processor has to treat the given fits file referenced by his number in the input list */
//int who_do_it(int size, int rank, int ii);

/*! this function determines whether the user list of detectors is correct or not */
void check_detector_is_in_fits(struct detectors det,struct detectors bolo_fits, std::string filename);

/*! this function determines whether reference and offsets HDUs are present and have the correct size */
int check_positionHDU(std::string fname,long ns,struct detectors det, int format, struct checkHDU &check_it);

/*! this function determines whether channels, time, signal and mask HDUs are present and have the correct size */
int check_commonHDU(std::string fname,long ns,struct detectors det, struct checkHDU &check_it);

/*! check RA/DEc table presence : only for HIPE format */
int check_altpositionHDU(std::string fname,long ns,struct detectors det, struct checkHDU &check_it);

/*! Check presence of non-flagged NANs in position tables */
int check_NAN_positionHDU(std::string fname,long ns,struct detectors det, struct checkHDU check_it);

/*! check presence of non-flagged NANs in time, signal and mask tables */
int check_NAN_commonHDU(std::string fname,long ns,struct detectors det, struct checkHDU check_it);

/*! check non-flagged NANs in RA/DEC HIPE format */
int check_NAN_altpositionHDU(std::string fname,long ns,struct detectors det, struct checkHDU check_it);

/*! Check that a given fits file bolometer list is exactly the same as the first input fits file one */
bool check_bolos(std::vector<std::string> bolo_fits_vect, std::vector<std::string> bolo_fits_0_vect);

/*!  Lookfor fully or more than 80% flagged detectors, also flag singletons */
int check_flag(std::string fname,struct detectors det,long ns, std::string outname,long *&bolos_global,long *&bolos_global_80, double *percent_tab, struct checkHDU check_it);

/*! check for time gaps in time table */
int check_time_gaps(std::string fname,long ns, double fsamp, struct common dir, struct checkHDU check_it);

/*! generating log files for user information */
void log_gen(long  *bolo_, std::string outname, struct detectors det, double *percent_tab=NULL);





#endif /* TOOLS_H_ */
