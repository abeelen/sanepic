#ifndef FITPOLY_H_
#define FITPOLY_H_

//! Finds the coefficients of a polynomial a(x) of degree n that fits the data sy : a(sx(i)) = sy(i)
/*!
 * WARNING  : taille must be > norder \n
 * sx is normalized between -1 and 1
 * \param norder The polynomial degree
 * \param taille sx and sy size is taille
 * \param sx A vector of "taille" points
 * \param sy Data to be fitted by "a"
 * \param a Row vector of length n+1 containing the polynomial coefficients in ascending powers
 */
void fitpoly(int norder, long taille, double *sx, double *sy, double *a);

#endif /* FITPOLY_H_ */
