#ifndef COORD_H_
#define COORD_H_

int convert_Dirfile_LON_LAT(struct samples samples_struct, struct param_sanePos pos_param, long iframe_min, long iframe_max, std::vector<std::vector<std::string> > bolo_vect);

void equatorial2galactic (long nx, double *RA, double *DEC, double **glon, double **glat);
void galactic2equatorial (long nx, double *glon, double *glat, double **RA, double **DEC);

void Spherical2Cartesian (double RA, double DEC, double Cart[3]);
void Cartesian2Spherical (double Cart[3], double *RA, double *DEC);

void RotVect (double Mat[3][3], double IN[3], double OUT[3]);
double range (double angle);
double ranrm (double angle);

#endif /* COORD_H_ */
