#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#include "localmathlib.h"

// ***************************************************************
// ********************** cDataInterpol Class
// ***************************************************************
dataInterpol::dataInterpol() { FreeData(); }

dataInterpol::~dataInterpol(void) { clear(); }

void dataInterpol::FreeData(void) {
	xval.clear();
	yval.clear();
	a0.clear();
	b0.clear();
	isempty_val = true;
}

void dataInterpol::clear()
{
	if (xval.size() > 0) xval.clear();
	if (yval.size() > 0) yval.clear();
	if (a0.size() > 0) a0.clear();
	if (b0.size() > 0) b0.clear();
	isempty_val = true;
}

void dataInterpol::set_points(const dvector &x,const dvector &y) {
	assert(x.size() == y.size());
	assert(x.size()>2);
	int n = (int)x.size();
	
	isempty_val = false;
	xval = x;
	yval = y;
	SortAscend(); // sort xval and yval

	xmin = xval.front(); ymin = yval.front();
	xmax = xval.back();  ymax = yval.back();

	double x1, x2, y1, y2;
	for (int i = 0; i < n - 1; i++) {
		a0.push_back(yval.at(i));
		x1 = xval.at(i);   y1 = yval.at(i);
		x2 = xval.at(i+1); y2 = yval.at(i+1);
		b0.push_back((y2- y1) / (x2 - x1));
	}
}

double dataInterpol::operator()(double x) {
	int ilo, n = (int)xval.size() - 2;

	// if x exceed the iteration domain [xmin, xmax] return end points [ymin, ymax] 
	if (x <= xmin) {
		return ymin;
	}
	else if (x >= xmax) {
		return ymax;
	}

	// find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
	dvector::iterator it;
	it = std::lower_bound(xval.begin(), xval.end(), x);
	ilo = (int)std::distance(xval.begin(), it) - 1;	    // pass lower bound index
	ilo = ilo > 0 ? ilo : 0;						// if ilo < 0 => ilo = 0
	ilo = ilo < n ? ilo : n - 1;					// if ilo > n => ilo = n
	double dx = x - xval.at(ilo);
	return a0.at(ilo) + b0.at(ilo)* dx;				//a + b*dx
}

void dataInterpol::SortAscend(void)
{
	int Ndata = (int)xval.size();
	sort_xy(&xval[0], &yval[0], sizeof(double), Ndata - 1);
}

bool dataInterpol::isempty(void) const {
	return isempty_val;
}


/***************************************************************************************
*******************              other miscelaneous functions                 ***********
****************************************************************************************/
double* linspace(double a, double b, unsigned int n) {
	double *output;

	output = new double[n];
	if (n < 2) {
		output[0] =  a;
		return output;
	}

	double step = (b - a) / (n - 1);
	for (unsigned int i = 0; i < n; i++){
		output[i] = a;
		a += step;
	}
	return output;
}

double* logspace(double a, double b, unsigned int n) {
	double *output;
	a = log10(a);
	b = log10(b);

	output = new double[n];
	if (n < 2) {
		output[0] = pow(10,a);
		return output;
	}

	double step = (b - a) / (n - 1);
	for (unsigned int i = 0; i < n; i++) {
		output[i] = pow(10, a);
		a += step;
	}
	return output;
}


int indexMin(double array[], int size)
{
	assert(array != NULL);
	assert(size >= 0);

	if ((array == NULL) || (size <= 0))
		return (int)NAN;

	double val = DBL_MAX;
	int iMin = -1;
	for (int i = 0; i < size; i++)
		if (array[i] < val) {
			iMin = i;
			val = array[i]; }
	return iMin;
}


double Min(double array[], int size)
{
	assert(array != NULL);
	assert(size >= 0);

	if ((array == NULL) || (size <= 0))
		return NAN;

	double val = -DBL_MAX;
	for (int i = 0; i < size; i++)
		if (val < array[i]) {
			val = array[i];
		}
	return val;
}

void sort_xy(double* x, void* y, const int y_size, int x_end, int x_front, bool ascend)
{
	int x_size = sizeof(double);
	// x_front = subscript of beginning of x array
	// x_end = subscript of end of x array
	if (x_front < 0) return;

	int middle;
	if (x_front < x_end)
	{
		middle = partition(x, y, y_size, x_end, x_front, ascend);
		sort_xy(x, y, y_size, middle - 1,    x_front, ascend);   // sort first section
		sort_xy(x, y, y_size,      x_end, middle + 1, ascend);   // sort second section
	}
	return;
}

int partition(double *x, void* y, const int y_size, int x_end, int x_front, bool ascend)
{
	int x_size = sizeof(double);
	double x0 = x[x_front];
	int i = x_front;
	int j = x_end + 1;
	do{
		do j--; while ((ascend ? x0 < x[j] : x0 > x[j]) && j > x_front);
		do i++; while ((ascend ? x0 > x[i] : x0 < x[i]) && i < x_end + 1);

		if (i < j) {
			swap((char *)x + i * x_size, (char *)x + j * x_size, x_size);
			swap((char *)y + i * y_size, (char *)y + j * y_size, y_size);
		}
	} while (i < j);
	swap((char *)x + (i-1)*x_size, (char *)x + x_front*x_size, x_size);
	swap((char *)y + (i-1)*y_size, (char *)y + x_front*y_size, y_size);
	
	return i - 1;           // returns middle subscript
}

void swap(void* vp1, void* vp2, const int size) {
	char *buffer = (char *)malloc(sizeof(char)*size);
	memcpy(buffer, vp1, size);
	memcpy(vp1, vp2, size);
	memcpy(vp2, buffer, size);
	free(buffer);
}

std::vector<std::string> Tokenize(std::string str, std::string braket){
	int i = 0;
	int w_end = 0;
	int w_begin = 0;
	bool openbraket = false;
	std::vector<std::string> tokens;
	const char *lbra = "(", *rbra = ")";

	//std::cout << "Tokenize function" << std::endl;
	//std::cout << "Text: " << str << std::endl;

	if (!braket.compare("(") || !braket.compare(")")) { lbra = "("; rbra = ")"; }
	if (!braket.compare("[") || !braket.compare("]")) { lbra = "["; rbra = "]"; }
	if (!braket.compare("{") || !braket.compare("}")) { lbra = "{"; rbra = "}"; }

	if (str.empty()) {
		tokens.push_back("");
		return tokens;
	}

	int state = 0;
	while (i < str.size())	{
		if (*lbra == str.at(i) && !openbraket) openbraket = true;
		if (*rbra == str.at(i) &&  openbraket) openbraket = false;
		
		if (!openbraket && 
			(' ' == str.at(i) || '\t' == str.at(i) || ',' == str.at(i)))
				w_begin++;
			
		else {
			w_end++;
			if (i == str.size() - 1) {
				//std::cout << "Word at: " << w_begin << " " << w_end << std::endl;
				tokens.push_back(str.substr(w_begin, w_end));
			}

			else {
				if (!openbraket && 
					(' ' == str.at(i + 1) || '\t' == str.at(i + 1) || ',' == str.at(i + 1))) {

					//std::cout << "Word at: " << w_begin <<  " " << w_end << std::endl;

					tokens.push_back(str.substr(w_begin, w_end));
					w_begin += w_end;
					w_end = 0;
				}
			}
		}
		i++;
	}
	return tokens;
}

std::vector<std::string> EvalFunction(std::string &fName, std::string braket){
	int i = 0;
	int w_end = 0;
	int w_begin = 0;
	std::vector<std::string> tokens;
	std::string arg;
	std::string lbra, rbra;

	if (!braket.compare("(") || !braket.compare(")")) { lbra = "("; rbra = ")"; }
	if (!braket.compare("[") || !braket.compare("]")) { lbra = "["; rbra = "]"; }
	if (!braket.compare("{") || !braket.compare("}")) { lbra = "{"; rbra = "}"; }

	if (fName.find(lbra) > fName.size() || fName.find(rbra) > fName.size()) {
		tokens.clear();
		return tokens;
	}
	else {
		// Get new string with arguments (no brakets)
		arg = fName.substr(fName.find(lbra) + 1,
			fName.find(rbra) - fName.find(lbra) - 1);
		
		// Get function name
		fName = fName.substr(0, fName.find(lbra));

		// read arguments (tokenize)
		tokens = Tokenize(arg);
	}

	// Give error for empty argument
	if (!tokens.empty()) {
		for (int i = 0; i < tokens.size(); i++) {
			if (tokens.at(i).empty()) {
				ErrorMsg("Error! Missing argument in " + fName);
			}
		}
	}
	else {
		ErrorMsg("Error! Missing argument in " + fName);
	}
	return tokens;
}
// effective media based on Bruggerman's model
// fV: volume fraction of particles
// eps1: dielectric constant of particles
// eps_host: dielectric constant of host
cdouble eff_Bruggerman(const double &fV, const cdouble &eps1, const cdouble &eps_host) {
	double f1 = fV, f2 = 1 - fV;
	cdouble eps2 = eps_host;
	cdouble eps_eff;

	cdouble f1_eps1 = (3.0 * f1 - 1)*eps1;
	cdouble f2_eps2 = (3.0 * f2 - 1)*eps2;
	cdouble sqrt_eps1eps2 = sqrt((f1_eps1 + f2_eps2)*(f1_eps1 + f2_eps2) + 8.0*eps1*eps2);

	eps_eff = 0.25*(f1_eps1 + f2_eps2 - sqrt_eps1eps2);
	
	if (eps_eff.imag() < 0 || (eps_eff.imag() < 1E-10 && eps_eff.real() < 0)) 
		eps_eff += 0.5*sqrt_eps1eps2;
	return eps_eff;
}

void ErrorMsg(std::string strError) {
	Log("Error! %s", strError.c_str());
	exit(EXIT_FAILURE);
}
 
/***************************************************************/
/* Simple general-purpose status logging (not thread-safe) *****/
/***************************************************************/
static std::string LogFileName = "";
static int LogToConsole = 0;

char* CurrentDate(void) {
	time_t curr_time;
	char *date;
	curr_time = time(NULL);
	date = ctime(&curr_time);
	date[strlen(date) - 1] = '\0'; // remove '\n' character at the end
	return date;
}

void SetConsoleLogging(void)
{
	LogToConsole = 1;
}

void SetLogFileName(std::string strlogFileName)
{
	LogToConsole = 0;
	if (!LogFileName.empty()) LogFileName = "";
	if (strlogFileName.empty()) strlogFileName = "logFile";
	strlogFileName += ".log";
	LogFileName = strlogFileName;
}

void Log(char const *format, ...)
{
	FILE *f = 0;

	// Construct formated message
	va_list ap;
	char buffer[VARARG_BUFLEN];
	va_start(ap, format);
	vsnprintf(buffer, VARARG_BUFLEN, format, ap);
	va_end(ap);

	// Print message
	if (LogToConsole)
		printf("%s. %s\n", CurrentDate(), buffer);

	else {
		f = fopen(LogFileName.c_str(), "a");
		if (f == 0) ErrorMsg("Unavailable log file");
		fprintf(f, "%s. %s\n", CurrentDate(), buffer);
		fclose(f);
	}
}