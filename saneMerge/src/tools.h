

#ifndef TOOLS_H_
#define TOOLS_H_

/*! this function tests that the detectors lists are the same in whole fits files, also tests that the files have a crescent time reference */
void file_compatibility_verification(struct samples samples_struct);

/*! copy reference position tables from each file to output file */
void copy_ref_pos(fitsfile *outfptr, struct samples samples_struct, long ns_final);

/*! copy offsets table from this file to output file */
void copy_offsets(fitsfile * fptr, fitsfile *outfptr);

/*! copy channels list from this file to output file */
void copy_channels(fitsfile * fptr, fitsfile *outfptr);

/*! copy time tables from each file to output file */
void copy_time(fitsfile *outfptr, struct samples samples_struct, long ns_final);

/*! copy signal tables from each file to output file */
void copy_signal(fitsfile *outfptr, struct samples samples_struct, struct detectors det, long ns_final);

/*! copy flag tables from each file to output file */
void copy_mask(fitsfile *outfptr, struct samples samples_struct, struct detectors det, long ns_final);

/*! copy RA and DEC tables (HIPE format only) from each file to output file */
void copy_RA_DEC(fitsfile *outfptr, struct samples samples_struct, struct detectors det, long ns_final);

#endif /* TOOLS_H_ */
