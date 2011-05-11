#ifndef CHOLESKY_H_
#define CHOLESKY_H_

//! Computes Cholesky factorisation for matrix "a"
/*!
 * a is symmetric, positive definite
 * l is a lower triangular matrix with strictly positive diagonal entries
 * \param n is the number of lines and columns in a
 * \param a is a squared matrix which side is n
 * \param l is the cholesky decomposition matrix that verify a = l * transpose(l)
 */
void cholesky(long n, double **a, double **l);

//! Resolve a triangular system L.y = b
/*!
 * We suppose that the system has a single solution
 * \param L is a n*m triangular matrix
 * \param n is the number of lines in L
 * \param m is the number of columns in L
 * \param y is the unknown array that has to be computed
 * \param b is an array of double
 * \param sup_inf A boolean that indicates whether L is upper triangular (sup_inf = 0) or lower triangular (sup_inf = 1)
 */
void system_triang(double **L, long n, long m, double *y, double *b, bool sup_inf);

//! Resolves the system a.x = b
/*!
 *  a : matrix symetric n*n, positive-definite
 * \param a is a squared matrix which side is n
 * \param b is an array of double that verify A.x = b
 * \param l is the matrix computed by cholesky routine
 * \param x is the unknown array that has to be computed
 * \param n is the number of lines and columns in a
 */
void solve_cholesky(double **a, double *b, double **l, double *x, long n);

#endif /* CHOLESKY_H_ */
