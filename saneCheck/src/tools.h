

#ifndef TOOLS_H_

void check_hdu(std::string fname,long ns, struct detectors det);
void check_NaN(std::string fname,long ns,struct detectors det);
int check_flag(std::string fvect, struct detectors det, long ns, std::string outname);

#define TOOLS_H_


#endif /* TOOLS_H_ */