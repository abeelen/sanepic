#ifndef TOOLS_H_
#define TOOLS_H_


/*! this function determines whether the user list of detectors is correct or not */
void check_detector_is_in_fits(std::vector<std::string> det, long ndet, std::vector<std::string> bolo_fits, std::string filename);

/*! this function determines whether reference and offsets HDUs are present and have the correct size */
int check_positionHDU(std::string fname,long ns, long ndet, int format, struct checkHDU &check_it);

/*! this function determines whether channels, time, signal and mask HDUs are present and have the correct size */
int check_commonHDU(std::string fname,long ns, long ndet, struct checkHDU &check_it);

/*! check RA/DEc table presence : only for HIPE format */
int check_altpositionHDU(std::string fname,long ns, long ndet, struct checkHDU &check_it);

/*! Check presence of non-flagged NANs in position tables */
int check_NAN_positionHDU(std::string fname,long ns, std::vector<std::string> det, long ndet, struct checkHDU check_it);

/*! check presence of non-flagged NANs in time, signal and mask tables */
int check_NAN_commonHDU(std::string fname,long ns, std::vector<std::string> det, long ndet, struct checkHDU check_it);

/*! check non-flagged NANs in RA/DEC HIPE format */
int check_NAN_altpositionHDU(std::string fname,long ns, std::vector<std::string> det, long ndet, struct checkHDU check_it);

/*! Check that a given fits file bolometer list is exactly the same as the first input fits file one */
bool check_bolos(std::vector<std::string> bolo_fits_vect, std::vector<std::string> bolo_fits_0_vect);

/*!  Lookfor fully or more than 80% flagged detectors, also flag singletons */
int check_flag(std::string fname, std::vector<std::string> det, long ndet, long ns, std::string outname, long *&bolos_global,
		long *&bolos_global_80, double *percent_tab, struct checkHDU check_it);

/*! check for time gaps in time table */
int check_time_gaps(std::string fname,long ns, double fsamp, struct param_common dir, struct checkHDU check_it);

/*! check whether the detector gain is correct or not and estimates the gain for each detector */
int check_bolo_gain(std::string fname,long ns, std::string bolo_gain_filename, std::vector<std::string> det, long ndet, struct checkHDU check_it);

/*! compute the median of a double vector */
double median(std::vector<double> vec);

/*! generating log files for user information */
void log_gen(long *bolo_, std::string outname, std::vector<std::string> det, long ndet, double *percent_tab=NULL);

/*! read the signal for all detector and for one given time sample */
int read_sample_signal_from_fits(string filename, int sample, double *& signal_samp, std::vector<std::string> det, long ndet);



#endif /* TOOLS_H_ */
