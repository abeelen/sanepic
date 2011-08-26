#ifndef TOOLS_H_
#define TOOLS_H_


/*! read saneCheck log files : get sample indices => where the gaps are */
int read_indices_file(std::string fname, struct param_common dir, std::vector<long> &indice, double &fsamp, long &init_num_delete, long &end_num_delete);

void refresh_indice(double &fsamp, long init_num_delete, long end_num_delete, std::vector <long> &indice, long ns);

/*!  compute the number of sample that must be added to fill the gaps and have a continous timeline  */
long how_many(std::string fname,long ns, std::vector <long> &indice, double fsamp, std::vector <long> &add_sample, std::vector <long> & suppress_time_sample);

/*! fix time table : draw a continuous time table by filling gaps using the calculated sampling frequency */
void fix_time(double *time, double *&time_fixed, std::vector <long> indice, std::vector <long> &add_sample, double fsamp, long nsamples_total, std::vector <long> suppress_time_sample, long init_num_delete);

/*! fix a row vector : given a row, fill this row with average values */
void fix_row(double *row, double *&row_fixed, std::vector <long> indice, std::vector <long> add_sample, long nsamples_total, std::vector <long> suppress_time_sample, long init_num_delete);

/*! Fill a mask row vector with ones for each gap generated values */
void fix_mask(int *mask, int *&mask_fixed, std::vector <long> indice, std::vector <long> add_sample, long nsamples_total, long init_num_delete);

/*! Copy time header from original fits file to output, insert final time table and update header */
void insert_time(fitsfile * fptr, fitsfile *outfptr, double *time, long ns_final);

/*! Insert one fixed row in the output fits file */
void insert_row_in_image(fitsfile *fptr, fitsfile *outfptr, std::string field, double *RA_fixed, long ns_total);

/*! insert a mask row in the output mask image */
void insert_mask_in_image(fitsfile *fptr, fitsfile *outfptr, std::string field, int *mask_fixed, long ns_total);

/*! insert "reference position" table in the output fits file */
void insert_ref_pos_in_fits(fitsfile *fptr, fitsfile *outfptr, double *RA, double *DEC,double *PHI, long ns_total);

/*! Copy input signal header to output and fill the gaps in signal table */
void fix_signal(fitsfile * fptr, fitsfile *outfptr, std::string name, long ns_total, std::vector<std::string> det, long ndet, std::vector <long> indice, std::vector<long> add_sample, std::vector <long> suppress_time_sample, long init_num_delete);

/*! Copy input RA and DEC header to output and fill the gaps in those tables : HIPE format only */
void fix_LON_LAT(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, std::vector<std::string> det, long ndet, std::vector <long> indice, std::vector<long> add_sample, std::vector <long> suppress_time_sample, long init_num_delete);

/*! Copy input mask header to output and fill the gaps with ones */
void fix_mask(fitsfile * fptr, fitsfile *outfptr, std::string name, long ns_total, std::vector<std::string> det, long ndet, std::vector <long> indice, std::vector<long> add_sample, std::vector <long> suppress_time_sample, long init_num_delete);

/*! Copy input time header to output and fill the gaps with computed values using sampling frequency */
void fix_time_table(fitsfile * fptr, fitsfile *outfptr, std::string name, long ns_total, std::vector <long> indice, std::vector<long> add_sample, long ns_origin, double fsamp, std::vector <long> suppress_time_sample, long init_num_delete);

/*! copy "reference position" table to output and fill the gaps with average values */
void fix_ref_pos(fitsfile * fptr, fitsfile *outfptr, std::string name, long ns_total, std::vector <long> indice, std::vector<long> add_sample, std::vector <long> suppress_time_sample, long init_num_delete);

#endif /* TOOLS_H_ */
