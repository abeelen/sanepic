#ifndef TOOLS_H_
#define TOOLS_H_

#include "struct_definition.h"

/*! Copy RA DEC and PHI (reference detector) tables from input fits to output */
void copy_ref_pos(fitsfile * fptr, fitsfile *outfptr, std::string name, long min_sample, long max_sample);

/*! copy bolometers offsets (in regards to ref detector) from input to output */
void copy_offsets(fitsfile * fptr, fitsfile *outfptr);

/*! Copy channel list table from input fits to output */
void copy_channels(fitsfile * fptr, fitsfile *outfptr);

/*! Copy resized time table from input fits to output */
void copy_time(fitsfile * fptr, fitsfile *outfptr, double *time, long min_sample, long max_sample);

/*! Copy resized signal table from input fits to output */
void copy_signal(fitsfile * fptr, fitsfile *outfptr,   std::string name, long min_sample, long max_sample, struct detectors det);

/*! Copy resized mask table from input fits to output */
void copy_mask(fitsfile * fptr, fitsfile *outfptr,   std::string name, long min_sample, long max_sample, struct detectors det);

/*! Copy resized RA and DEC tables from input fits to output */
void copy_RA_DEC(fitsfile * fptr, fitsfile *outfptr, string name, long min_sample, long max_sample, struct detectors det);

#endif /* TOOLS_H_ */
