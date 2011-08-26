#ifndef TOOLS_H_
#define TOOLS_H_

#include "struct_definition.h"

//! Copy reference position table from input fits to output fits file
/*!
 * Copy only needed samples : contained in [min_sample, max_sample]
 \param fptr A pointer to the input fits file
 \param outfptr A pointer to the output fits file in which the reference position table will be copied
 \param name Input fits filenmae
 \param min_sample First sample index : copy samples from min_sample
 \param max_sample Last sample index : copy samples until max_sample
 */
void copy_ref_pos(fitsfile * fptr, fitsfile *outfptr, std::string name, long min_sample, long max_sample);

//! Copy resized time table from input fits to output fits file
/*!
 * Copy only needed samples : contained in [min_sample, max_sample]
 \param fptr A pointer to the input fits file
 \param outfptr A pointer to the output fits file in which the time table will be copied
 \param time Input time array (read in input fits time table)
 \param min_sample First sample index : copy samples from min_sample
 \param max_sample Last sample index : copy samples until max_sample
 */
void copy_time(fitsfile * fptr, fitsfile *outfptr, double *time, long min_sample, long max_sample);

//! Copy resized signal table from input fits to output fits file
/*!
 * Copy only needed samples : contained in [min_sample, max_sample]
 * Copy table channel by channel
 \param fptr A pointer to the input fits file
 \param outfptr A pointer to the output fits file in which the signal table will be copied
 \param name Input fits filenmae
 \param min_sample First sample index : copy samples from min_sample
 \param max_sample Last sample index : copy samples until max_sample
 \param det A channel list, read from input fits file
 \param ndet det number of channels
 */
void copy_signal(fitsfile * fptr, fitsfile *outfptr, std::string name, long min_sample, long max_sample, std::vector<std::string> det, long ndet);

//! Copy resized mask table from input fits to output fits file
/*!
 * Copy only needed samples : contained in [min_sample, max_sample]
 * Copy table channel by channel
 \param fptr A pointer to the input fits file
 \param outfptr A pointer to the output fits file in which the mask table will be copied
 \param name Input fits filenmae
 \param min_sample First sample index : copy samples from min_sample
 \param max_sample Last sample index : copy samples until max_sample
 \param det A channel list, read from input fits file
 \param ndet det number of channels
 */
void copy_mask(fitsfile * fptr, fitsfile *outfptr, std::string name, long min_sample, long max_sample, std::vector<std::string> det, long ndet);

//! Copy resized LON and LAT tables from input fits to output fits file
/*!
 * Copy only needed samples : contained in [min_sample, max_sample]
 * Copy tables channel by channel
 \param fptr A pointer to the input fits file
 \param outfptr A pointer to the output fits file in which the RA/DEC tables will be copied
 \param name Input fits filenmae
 \param min_sample First sample index : copy samples from min_sample
 \param max_sample Last sample index : copy samples until max_sample
 \param det A channel list, read from input fits file
 \param ndet det number of channels
 */
void copy_LON_LAT(fitsfile * fptr, fitsfile *outfptr, string name, long min_sample, long max_sample, std::vector<std::string> det, long ndet);

#endif /* TOOLS_H_ */
