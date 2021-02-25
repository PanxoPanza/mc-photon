#pragma once
#include <assert.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <complex>
#include <time.h>

#include <cfloat>
#include <cstring>
#include <cstdio>
#include <stdarg.h>

#ifdef __linux__ 
//linux code goes here
#include "uni_localmathlib.h"
#elif _WIN32
// windows code goes here
#include "win_localmathlib.h"
#else

#endif

#define VARARG_BUFLEN 1000

typedef std::vector<double> dvector;
typedef std::complex<double> cdouble;

double* linspace(double a, double b, unsigned int n);
double* logspace(double a, double b, unsigned int n);
int indexMin(double array[], int size); //Return index of minimum
double Min(double array[], int size); //Return minimum value
void sort_xy(double* x, void* y, const int y_size, int x_end, int x_front = 0, bool ascend = true);
int partition(double* x, void* y, const int y_size, int x_end, int x_front = 0, bool ascend = true);
std::vector<std::string> Tokenize(std::string str, std::string braket = "(");
std::vector<std::string> EvalFunction(std::string &fName, std::string braket = "(");
void ErrorMsg(std::string strError);
cdouble eff_Bruggerman(const double &fV, const cdouble &eps1, const cdouble &eps_host);
void swap(void* vp1, void* vp2, const int size);
//vector<std::string> tokenize(const std::string &text);

class dataInterpol {
private:
	double xmin, xmax, ymin, ymax;
	dvector xval,  yval;
	dvector a0,  b0;
	bool isempty_val;

	void FreeData(void);
	void SortAscend(void);

public:
	dataInterpol();
	~dataInterpol();
	void set_points(const dvector &x, const dvector &y);
	double operator() (double w0);
	bool isempty(void) const;
	void clear();
};

/***************************************************************/
/* General-purpose status logging ******************************/
/***************************************************************/
void SetLogFileName(std::string logFileName);
void SetConsoleLogging(void);
void Log(const char *format, ...);
char* CurrentDate(void);