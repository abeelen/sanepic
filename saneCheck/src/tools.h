

#ifndef TOOLS_H_

int check_path(std::string strPath, std::string path_type);
int who_do_it(int size, int rank, int ii);
void check_detector_is_in_fits(struct detectors det,struct detectors bolo_fits, std::string filename);
void check_positionHDU(std::string fname,long ns,struct detectors det, int format);
void check_commonHDU(std::string fname,long ns,struct detectors det);
void check_altpositionHDU(std::string fname,long ns,struct detectors det);
void check_NAN_positionHDU(std::string fname,long ns,struct detectors det);
void check_NAN_commonHDU(std::string fname,long ns,struct detectors det);
void check_NAN_altpositionHDU(std::string fname,long ns,struct detectors det);
void check_NaN(std::string fname,long ns,struct detectors det);
bool check_bolos(std::vector<std::string> bolo_fits_vect, std::vector<std::string> bolo_fits_0_vect);
void check_flag(std::string fname,struct detectors det,long ns, std::string outname,long *&bolos_global,long *&bolos_global_80);
void check_time_gaps(std::string fname,long ns, double fsamp, struct common dir);
void log_gen(long  *bolo_, std::string outname, struct detectors det);


#define TOOLS_H_


#endif /* TOOLS_H_ */
