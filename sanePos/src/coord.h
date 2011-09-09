/*
 * coord.h
 *
 *  Created on: 9 sept. 2011
 *      Author: abeelen
 */

#ifndef COORD_H_
#define COORD_H_

void equatorial2galactic (long nx, double *RA, double *DEC, double **glon, double **glat);
void galactic2equatorial (long nx, double *glon, double *glat, double **RA, double **DEC);

void Spherical2Cartesian (double RA, double DEC, double Cart[3]);
void Cartesian2Spherical (double Cart[3], double *RA, double *DEC);

void RotVect (double Mat[3][3], double IN[3], double OUT[3]);
double range (double angle);
double ranrm (double angle);


#endif /* COORD_H_ */
