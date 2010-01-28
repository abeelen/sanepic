

#ifndef TOOLS_H_

void read_bolo_list(string fname, struct detectors &det);
void check_hdu(std::string fname,long ns, struct detectors det);
void check_NaN(std::string fname,long ns,struct detectors det);
int check_flag(std::string fvect, struct detectors det, long ns, std::string outname,std::vector<std::string> &bolos_global,std::vector<std::string> &bolos_global_80);
void check_time_gaps(string fname,long ns);
void log_gen(std::vector<std::string> &bolo_, std::string outname);


#define TOOLS_H_


#endif /* TOOLS_H_ */
