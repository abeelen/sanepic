#ifndef TOOLS_H_
#define TOOLS_H_


//! A structure that is used to sort a doubles vectors with standard routine "sort"
struct sortclass_double {
	bool operator() (double i,double j) { return (i<j);}
};

//! Determines whether the user list of detectors is correct or not
/*!
 * Prints a WARNING message to screen in case some differences are found
 \param det A channel list, read from bolovect lists (generated using ini file informations)
 \param ndet det's number of channels
 \param bolo_fits A channel list, read from fits file
 \param filename The considered scan file name (used to output warning message)
 */
void check_detector_is_in_fits(std::vector<std::string> det, long ndet, std::vector<std::string> bolo_fits, std::string filename);

//! Checks whether reference and offsets HDUs are present and have the correct size, and fills checkHDU structure in case of non-attendance
/*!
 * Prints a WARNING message to screen in case a table is missing. Return error code -1 in case tables are present but not conform (wrong column numbers for Reference positions for example)
 \param fname The considered scan file name
 \param ns This scans's number of samples
 \param ndet number of detectors in fits channels table (is used as a dimension for some tables)
 \param format An integer specifying "fname" file format : HIPE (1) or SANEPIC (2)
 \param check_it A checkHDU structure that determines which tables are missing
 \return error code -1 or 0 (everything OK)
 */
int check_positionHDU(std::string fname, long ns, long ndet, int format, struct checkHDU &check_it);

//! Determines whether channels, time, signal and mask HDUs are present and have the correct sizes, and fills checkHDU structure in case of non-attendance
/*!
 * Prints a WARNING message to screen in case a table is missing. Return error code -1 in case tables are present but not conform (wrong number of elements for time table for example)
 \param fname The considered scan file name
 \param ns This scans's number of samples
 \param ndet number of detectors in fits channels table (is used as a dimension for some tables)
 \param check_it A checkHDU structure that determines which tables are missing
 \return error code -1 or 0 (everything OK)
 */
int check_commonHDU(std::string fname, long ns, long ndet, struct checkHDU &check_it);

//! check LON/LAT table presence : only for HIPE format, and fills checkHDU structure in case of non-attendance
/*!
 * Prints a WARNING message to screen in case a table is missing. Return error code -1 in case tables are present but not conform (wrong number of elements for time table for example)
 \param fname The considered scan file name
 \param ns This scans's number of samples
 \param ndet number of detectors in fits channels table (is used as a dimension for some tables)
 \param check_it A checkHDU structure that determines which tables are missing
 \return error code -1 or 0 (everything OK)
 */
int check_altpositionHDU(std::string fname, long ns, long ndet, struct checkHDU &check_it);

//! Check presence of non-flagged NANs in position tables
/*!
 * Prints a WARNING message to screen in case a NAN has been found : Prints table name, and index in this table
 \param fname The considered scan file name
 \param ns This scans's number of samples
 \param det A channel list, read from fits file
 \param ndet number of detectors in fits channels table (is used as a dimension for some tables)
 \param check_it A checkHDU structure that determines which tables has to be checked according to check_*HDU routines
 \return The total number of NANs found
 */
long check_NAN_positionHDU(std::string fname, long ns, std::vector<std::string> det, long ndet, struct checkHDU check_it);

//! check presence of non-flagged NANs in time, signal and mask tables
/*!
 * Prints a WARNING message to screen in case a NAN has been found : Prints table name, and index in this table
 \param fname The considered scan file name
 \param ns This scans's number of samples
 \param det A channel list, read from fits file
 \param ndet number of detectors in fits channels table (is used as a dimension for some tables)
 \param check_it A checkHDU structure that determines which tables has to be checked according to check_*HDU routines
 \return The total number of NANs found
 */
long check_NAN_commonHDU(std::string fname, long ns, std::vector<std::string> det, long ndet, struct checkHDU check_it);

//! check persence of non-flagged NANs in LON/LAT HIPE format
/*!
 * Prints a WARNING message to screen in case a NAN has been found : Prints table name, and index in this table
 \param fname The considered scan file name
 \param ns This scans's number of samples
 \param det A channel list, read from fits file
 \param ndet number of detectors in fits channels table (is used as a dimension for some tables)
 \param check_it A checkHDU structure that determines which tables has to be checked according to check_*HDU routines
 \return The total number of NANs found
 */
long check_NAN_altpositionHDU(std::string fname, long ns, std::vector<std::string> det, long ndet, struct checkHDU check_it);

//! Check that a given fits file bolometer list is exactly the same as the first input fits file one
/*!
 * Prints a WARNING message to screen in case a NAN has been found : Prints table name, and index in this table
 \param bolo_fits_vect A channel list, read from current fits file
 \param bolo_fits_0_vect A channel list, read from first input fits file (according to samples_struct and fits_filelist.txt)
 \return A boolean : True if list are different, False otherwise
 */
bool check_bolos(std::vector<std::string> bolo_fits_vect, std::vector<std::string> bolo_fits_0_vect);

//!  Lookfor fully or more than 80% flagged detectors, also flag singletons.
/*!
 * Fills bolos_global, n_hund, bolos_global_80, n_heig, init_flag_num, end_flag_num and percent_tab
 \param fname The considered scan file name
 \param ns This scans's number of samples
 \param det A channel list, read from fits file
 \param ndet number of detectors in fits channels table (is used as a dimension for some tables)
 \param bolos_global A vector of binary flags, which sizes is ndet. if bolos_global[x] > 0 : det[x] is more than 99% flagged
 \param n_hund The number of detectors more than 99% flagged (also the number of bolos_global values >0)
 \param bolos_global_80 A vector of binary flags, which sizes is ndet. if bolos_global_80[x] > 0 : det[x] is more than 80% flagged
 \param n_heig The number of detectors more than 80% flagged (also the number of bolos_global_80 values >0)
 \param percent_tab An array containing the percentage of flagged samples for the channels (more than 80% only, to be printed in output log)
 \param init_flag_num If a number of samples ni is found at the beginning of whole flag table (for every detectors, in fits file), init_flag_num is set to ni
 \param end_flag_num If a number of samples ne is found at the end of whole flag table (for every detectors, in fits file), end_flag_num is set to ne
 \param check_it A checkHDU structure that determines which tables has to be checked according to check_*HDU routines
 \return An integer specifying if there were an error (>0) or not (=0)
 */
int check_flag(string fname, std::vector<std::string> det, long ndet, long ns,
		double *percent_tab, long &init_flag_num, long &end_flag_num);

//! Check for time gaps in fits time table
/*!
 * Locates time gaps and computes most populated frequency (to avoid dealing with gaps)
 * Fills indice and Populated_freq
 \param fname The considered scan file name
 \param ns This scans's number of samples
 \param fsamp Instrument Sampling frequency
 \param indice Sample indice : where a time gap has been found
 \param Populated_freq The most populated frequency only is kept and compare to the user frequency given in ini file (fsamp)
 \param check_it A checkHDU structure that determines which tables has to be checked according to check_*HDU routines
 \return An integer specifying if there were an error (>0) or not (=0)
 */
int check_time_gaps(std::string fname, long ns, double fsamp, std::vector<long> &indice, double &Populated_freq, struct checkHDU check_it);

/* check whether the detector gain is correct or not and estimates the gain for each detector */
// not used ...
int check_bolo_gain(std::string fname, long ns, std::string bolo_gain_filename, std::vector<std::string> det, long ndet, struct checkHDU check_it);

//! Computes the median of a double vector
/*!
 \param vec A vector of doubles
 \return A double corresponding to the median value of vec
 */
double median(std::vector<double> vec);

//!  Save informations to disk for saneFix program
/*!
 \param tmp_dir A string containing the temporary files pathname
 \param filename The considered scan file name
 \param init_flag Number of beginning flagged samples (present in each signal data), that has to be erased from fits tables by saneFix
 \param end_flag Number of ending flagged samples (present in each signal data), that has to be erased from fits tables by saneFix
 \param Populated_freq The most populated frequency (without time gaps)
 \param indice Time gaps indices : to be found and filled by saneFix
 \return An integer specifying if there were an error (>0) or not (=0)
 */
int print_to_bin_file(std::string tmp_dir, std::string filename, long init_flag, long end_flag, double Populated_freq, std::vector<long> indice);

//! generating log file (for user information), containing bad channels and their flagged percentage
/*!
 \param bolo_ A vector of binary flags, which sizes is ndet. Could be bolos_global_80 or bolos_global
 \param outname Log file output name
 \param det A channel list, read from fits file
 \param ndet number of detectors in fits channels table (is used as a dimension for some tables)
 \param percent_tab An array containing the percentage of flagged samples for the bad channels (more than 80% only, so it is an optional field)
 \return An integer specifying if there were an error (>0) or not (=0)
 */
void log_gen(long *bolo_, std::string outname, std::vector<std::string> det, long ndet, double *percent_tab=NULL);

/* read the signal for all detector and for one given time sample */
// Called by check_bolo_gain, not used ...
int read_sample_signal_from_fits(string filename, int sample, double *& signal_samp, std::vector<std::string> det, long ndet);



#endif /* TOOLS_H_ */
