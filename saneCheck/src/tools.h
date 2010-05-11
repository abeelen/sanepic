

#ifndef TOOLS_H_

/*! this function determines which processor has to treat the given fits file referenced by is number in the input list */
int who_do_it(int size, int rank, int ii);

/*! this function determines whether the user list of detectors is correct or not */
void check_detector_is_in_fits(struct detectors det,struct detectors bolo_fits, std::string filename);

/*! this function determines whether reference and offsets HDUs are present and have the correct size */
void check_positionHDU(std::string fname,long ns,struct detectors det, int format);

/*! this function determines whether channels, time, signal and mask HDUs are present and have the correct size */
void check_commonHDU(std::string fname,long ns,struct detectors det);

/*! check RA/DEc table presence : only for HIPE format */
void check_altpositionHDU(std::string fname,long ns,struct detectors det);

/*! Check presence of non-flagged NANs in position tables */
void check_NAN_positionHDU(std::string fname,long ns,struct detectors det);

/*! check presence of non-flagged NANs in time, signal and mask tables */
void check_NAN_commonHDU(std::string fname,long ns,struct detectors det);

/*! check non-flagged NANs in RA/DEC HIPE format */
void check_NAN_altpositionHDU(std::string fname,long ns,struct detectors det);

/*! Check that a given fits file bolometer list is exactly the same as the first input fits file one */
bool check_bolos(std::vector<std::string> bolo_fits_vect, std::vector<std::string> bolo_fits_0_vect);

/*!  Lookfor fully or more than 80% flagged detectors, also flag singletons */
void check_flag(std::string fname,struct detectors det,long ns, std::string outname,long *&bolos_global,long *&bolos_global_80);

/*! check for time gaps in time table */
void check_time_gaps(std::string fname,long ns, double fsamp, struct common dir);

/*! generating log files for user information */
void log_gen(long  *bolo_, std::string outname, struct detectors det);


#define TOOLS_H_


#endif /* TOOLS_H_ */
