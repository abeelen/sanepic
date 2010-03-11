#ifndef SANEINVIO_H_
#define SANEINVIO_H_

#include <string>
#include <vector>



/*!  This function writes the Inverse Covariance Matrices in binary format */
void write_InvNoisePowerSpectra(std::vector<std::string> bolos, long nbins, double * ell, double **Rellth, std::string ouputDir, std::string suffix);

/*! This function reads the NoiseNoise Matrices in a fits file (also reads the mixing matrices) */
void read_CovMatrix(std::string fname, std::vector<std::string> &bolos, long &nbins, double *&ell, double **&Rellth);


#endif /* SANEINVIO_H_ */
