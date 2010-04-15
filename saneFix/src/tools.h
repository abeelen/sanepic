#ifndef TOOLS_H_
#define TOOLS_H_

int read_indices_file(std::string fname, struct common dir, std::vector<long> &indice, double &fsamp);
//void read_tables_till_indice(std::string fname, long marker);
long how_many(std::string fname,long ns, std::vector <long> indice, double *&time, double fsamp, std::vector <long> &add_sample);
void fix_time(double *time, double *&time_fixed, std::vector <long> indice, std::vector <long> &add_sample, double fsamp, long nsamples_total);
void fix_row(double *row, double *&row_fixed, std::vector <long> indice, std::vector <long> add_sample, long nsamples_total);
void fix_mask(int *mask, int *&mask_fixed, std::vector <long> indice, std::vector <long> add_sample, long nsamples_total);
void insert_time(fitsfile * fptr, fitsfile *outfptr, double *time, long ns_final);
void insert_row_in_image(fitsfile *fptr, fitsfile *outfptr, std::string field, double *RA_fixed, long ns_total);
void insert_mask_in_image(fitsfile *fptr, fitsfile *outfptr, std::string field, int *mask_fixed, long ns_total);
void insert_ref_pos_in_fits(fitsfile *fptr, fitsfile *outfptr, double *RA, double *DEC,double *PHI, long ns_total);
void copy_offsets(fitsfile * fptr, fitsfile *outfptr);
void copy_channels(fitsfile * fptr, fitsfile *outfptr);
void fix_signal(fitsfile * fptr, fitsfile *outfptr, std::string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample );
//void fix_RA(fitsfile * fptr, fitsfile *outfptr, std::string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample );
//void fix_DEC(fitsfile * fptr, fitsfile *outfptr, std::string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample );
void fix_RA_DEC(fitsfile * fptr, fitsfile *outfptr, string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample);
void fix_mask(fitsfile * fptr, fitsfile *outfptr, std::string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample );
void fix_time_table(fitsfile * fptr, fitsfile *outfptr, std::string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample, long ns_origin, double fsamp);
void fix_ref_pos(fitsfile * fptr, fitsfile *outfptr, std::string name, long ns_total, struct detectors det, std::vector <long> indice, std::vector<long> add_sample );
#endif /* TOOLS_H_ */
