#include <complex>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include "localmathlib.h"

#define SPEEDOFLIGHT 299792458 // m/s
#define ECHARGE 1.602176565e-19 // C
#define EPS0 8.85418782e-12 // F/m
#define MU0 M_PI*4e-7
#define H_PLANK 6.62606957e-34
#define hbar H_PLANK / (2.*M_PI)

using namespace std;

typedef complex<double> cdouble;
typedef vector<double> dvector;
typedef vector<cdouble> cvector;
typedef cdouble(*epsmu_fn)(double);
typedef vector<string> vstring;

class cDataInterpol {
private:
	dvector xdata, zRe, zIm;
	cvector a0, b0;
	void CleanUp(void);
	void Initiallize(void);
	void xBubbleSort(void);

public:
	double xmin, xmax;
	cdouble zmin, zmax;
	cDataInterpol();
	cDataInterpol(string File, bool is_application_path = true);
	~cDataInterpol();
	string get_path(bool is_application_path);
	bool GetFromFile(string File, bool is_application_path = true);
	cdouble operator() (double w0);
	};

class Material
{
private:
	cDataInterpol eps_interp, mu_interp;
	epsmu_fn eps, mu;
	cdouble eps0;
	string MatType;
	//vstring EvalFunction(string &fName);

public:
	Material();
	Material(const string &MaterialSet);
	~Material();
	int SetMaterial(const string &MaterialSet);
    string MaterialType(void) const;
	cdouble Eps(double w);
	cdouble Mu(double w);
};

void printEps(epsmu_fn fx, double w1, double w2, int Nw);
cdouble epsAu(double w);
cdouble epsAg(double w);
cdouble epsAl(double w);
cdouble epsAlN(double w);
cdouble epsCr(double w);
cdouble epsW(double w);
cdouble epsPt(double w);
cdouble epsSiC(double w);
cdouble epsSiDoped(double w, double N=0.);
cdouble epsSiO2(double w);
cdouble epsVO2h(double w);
cdouble epsVO2c(double w);
cdouble epsMoO3(double w);
cdouble epsKBr(double w);
cdouble epsAl2O3(double w);
cdouble epsSi3N4(double w);
cdouble epsFePt(double w);
cdouble epsC(double w);
cdouble epsBN(double w); // c-BN
cdouble epsVac(double w);
cdouble muVac(double w);
cdouble epsPDMS(double w);
cdouble epsPMMA(double w);
cdouble epsHDPE(double w);
cdouble epsLDPE(double w);
cdouble epsPVdF_HFP(double w);
cdouble epsConst(double w, cdouble eps);
cdouble epsDrude(double w, double wp, double gamma, double epsInf = 1);
//#endif
