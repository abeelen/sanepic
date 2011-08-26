#ifndef TOOLS_H_
#define TOOLS_H_

//! Test that the detectors lists are the same in whole fits files, also test that the files have a crescent time reference
/*!
 * No return Value, Exit is called in case of incompatibility
 \param samples_struct A structure that contains input fits filenames
 */
void file_compatibility_verification(string dirfile, struct samples samples_struct);

//! copy reference position table from each file to output file
/*!
 \param outfptr A pointer to the output opened fits file in which the reference position table will be copied
 \param samples_struct A structure that contains input fits filenames
 \param ns_final ouput fits file number of samples (sum of input number of samples)
 */
void copy_ref_pos(fitsfile *outfptr, string dirfile, struct samples samples_struct, long ns_final);

//! copy time table from each file to output file
/*!
 \param outfptr A pointer to the output opened fits file in which the time table will be copied
 \param samples_struct A structure that contains input fits filenames
 \param ns_final ouput fits file number of samples (sum of input number of samples)
 */
void copy_time(fitsfile *outfptr, string dirfile, struct samples samples_struct, long ns_final);

//! copy signal table from each file to output file
/*!
 * Copy tables detector by detector
 \param outfptr A pointer to the output opened fits file in which the signal table will be copied
 \param samples_struct A structure that contains input fits filenames
 \param det A channel list, read from input fits file
 \param ndet det number of channels
 \param ns_final ouput fits file number of samples (sum of input number of samples)
 */
void copy_signal(fitsfile *outfptr, string dirfile, struct samples samples_struct, std::vector<std::string> det, long ndet, long ns_final);

//! copy flag table from each file to output file
/*!
 * Copy tables detector by detector
 \param outfptr A pointer to the output opened fits file in which the signal table will be copied
 \param samples_struct A structure that contains input fits filenames
 \param det A channel list, read from input fits file
 \param ndet det number of channels
 \param ns_final ouput fits file number of samples (sum of input number of samples)
 */
void copy_mask(fitsfile *outfptr, string dirfile, struct samples samples_struct, std::vector<std::string> det, long ndet, long ns_final);

//! copy RA and DEC tables (HIPE format only) from each file to output file
/*!
 * Copy tables detector by detector
 \param outfptr A pointer to the output opened fits file in which the signal table will be copied
 \param samples_struct A structure that contains input fits filenames
 \param det A channel list, read from input fits file
 \param ndet det number of channels
 \param ns_final ouput fits file number of samples (sum of input number of samples)
 */
void copy_LON_LAT(fitsfile *outfptr, string dirfile, struct samples samples_struct, std::vector<std::string> det, long ndet, long ns_final);

#endif /* TOOLS_H_ */
