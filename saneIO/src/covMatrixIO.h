//#ifndef COVMATRIX_H_
//#define COVMATRIX_H_
//
//
//#include <string>
//#include <vector>
//



//void read_noisefile(string fname, string bolo1bolo2, double *ell, double *SPN, long *nbins);

///*! Return the longest string length of a string vector */
//long maxStringLength(std::vector<std::string> strings);


///*!  This function writes the Inverse Covariance Matrices in binary format */
//void write_InvNoisePowerSpectra(std::vector<std::string> bolos, long nbins, double * ell, double **Rellth, std::string ouputDir, std::string suffix);

///*! This function reads the Inverse Covariance Matrices in binary format */
//void read_InvNoisePowerSpectra(std::string prefix, std::string boloName, std::string suffix, long * nbins, long * ndet, double ** ell, double *** SpN_all);

///*! This function writes the NoiseNoise Matrices in a fits file (also Writes the mixing matrices) */
//void write_CovMatrix(std::string fname, std::vector<std::string> bolos, long nbins, double *ell, double **Rellth);
//void write_CovMatrix2(std::string fname, std::vector<std::string> bolos, long nbins, double *ell, double **Rellth);

/*! Writes the reduced mixing matrix in a binary file */
//void write_ReducedMixingMatrix(double **mixmat,long ndet,int ncomp, std::string outputDir);

///*! Writes the reduced mixing matrix from a binary file */
//void read_ReducedMixingMatrix(double **&mixmat,long &ndet,int &ncomp, std::string dir);
//
///*! This function reads the NoiseNoise Matrices in a fits file (also reads the mixing matrices) */
//void read_CovMatrix(std::string fname, std::vector<std::string> &bolos, long &nbins, double *&ell, double **&Rellth);
////
/////*! Returns the table Format for the given vector of string*/
////char* tableFormat(std::vector<std::string> strings);
////
/////*! Transform a vector of string into a array of char*/
////char** vString2carray(std::vector<std::string> strings);
//
//#endif /* COVMATRIX_H_ */
