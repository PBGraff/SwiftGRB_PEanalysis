#ifndef __BAMBI_H__
#define __BAMBI_H__ 1

void FirstRunCheck(int);
void GetNetTol();
void CopyFile(std::string, std::string);

void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void bambi(int &ndata, int &ndim, double **BAMBIData, double &lowlike);

void getLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void getphysparams(double *Cube, int &ndim, int &nPar, void *context);
void getallparams(double *Cube, int &ndim, int &nPar, void *context);

#endif
