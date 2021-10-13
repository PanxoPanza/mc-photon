#define _USE_MATH_DEFINES
#include "materials.h"

#define OPTFOLDER "/optical_properties/"

typedef complex<double> cdouble;
typedef cdouble (*epsmu_fn)(double);
const cdouble II(0.0,1.0);

// ***************************************************************
// ********************** cDataInterpol Class
// ***************************************************************
cDataInterpol::cDataInterpol() { Initiallize(); }

cDataInterpol::cDataInterpol(string File, bool is_application_path)
{
	Initiallize();
	GetFromFile(File, is_application_path);
}

cDataInterpol::~cDataInterpol(void) { CleanUp(); }

void cDataInterpol::Initiallize(void) {
	xdata.clear();
	zRe.clear();
	zIm.clear();
	a0.clear();
	b0.clear();
}

void cDataInterpol::CleanUp(void)
{
	Initiallize();
}

string cDataInterpol::get_path(bool is_application_path){
		string s_path;
		if (is_application_path) // if the file is in the application directory
			s_path = get_application_path() + OPTFOLDER;

		else
			s_path = get_working_path(); // if file is in working directory

		return s_path;
}

bool cDataInterpol::GetFromFile(string File, bool is_application_path) {
	// Get full bath + FileName
	File = get_path(is_application_path) + File;

	ifstream infile(File);
	// Read data from file
	if (infile.is_open()) { // only if File can be opened
		double data;
		int Ndata;

		// Extract data (cross sections in mm^2)
		infile >> Ndata;

		//while(!infile.eof()){
		for (int i = 0; i < Ndata; i++) {
			if (infile.eof())
				ErrorMsg("Error! Reach end of file before completing reading of data");

			infile >> data; xdata.push_back(data);
			infile >> data; zRe.push_back(data);
			infile >> data; zIm.push_back(data);
			//cout << xdata.at(i) << ' ' << zRe.at(i) << ' ' << zIm.at(i) << endl;
		}
		infile.close();

		xBubbleSort(); // sort data in ascending order

		// set interpolation parameters
		double x1, x2;
		cdouble z1, z2;
		xmin = xdata.front(); zmin = zRe.front() + II*zIm.front();
		xmax = xdata.back();  zmax = zRe.back()  + II*zIm.back();
		for (int i = 0; i < Ndata - 1; i++) {
			a0.push_back(zRe.at(i) + II*zIm.at(i));
			x1 = xdata.at(i);    z1 = zRe.at(i)   + II*zIm.at(i);
			x2 = xdata.at(i+1);  z2 = zRe.at(i+1) + II*zIm.at(i+1);
			b0.push_back((z2 - z1) / (x2 - x1));
		}
	}
	else {
		ErrorMsg("Error! Unable to open file: " + File);
	}
	return true;
}

cdouble cDataInterpol::operator()(double x) {

	if		(x <= xmin) return zmin;
	else if (x >= xmax) return zmax;
	else {
		// find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
		dvector::iterator it;
		it = lower_bound(xdata.begin(), xdata.end(), x);
		int ilo = max(int(it - xdata.begin()) - 1, 0);
		double dx = x - xdata.at(ilo);
		return a0.at(ilo) + b0.at(ilo)* dx;
	}
}

void cDataInterpol::xBubbleSort(void)
{
	int i, j, flag = 1;    // set flag to 1 to start first pass
	double temp;             // holding variable
	int numLength = (int)xdata.size();
	for (i = 1; (i <= numLength) && flag; i++)
	{
		flag = 0;
		for (j = 0; j < (numLength - 1); j++)
		{
			if (xdata.at(j) > xdata.at(j + 1))      // ascending order
			{
				temp = xdata.at(j);             // swap elements
				xdata.at(j) = xdata.at(j + 1);
				xdata.at(j + 1) = temp;

				temp = zRe.at(j);             // swap elements
				zRe.at(j) = zRe.at(j + 1);
				zRe.at(j + 1) = temp;

				temp = zIm.at(j);             // swap elements
				zIm.at(j) = zIm.at(j+1);
				zIm.at(j + 1) = temp;
				flag = 1;               // indicates that a swap occurred.
			}
		}
	}
	return;   //arrays are passed to functions by address; nothing is returned
}

// end of cDataInterpol Class

// ***************************************************************
// ************************* Material Class
// ***************************************************************
Material::Material()
{
}

Material::Material(const string &MaterialSet)
{
	eps = NULL;
	mu = NULL;
	SetMaterial(MaterialSet);
}

Material::~Material()
{
	if (eps != NULL) eps = NULL;
	if (mu != NULL) mu = NULL;
}

int Material::SetMaterial(const string &MaterialSet)
{
	MatType = MaterialSet;
	vstring Tokens = EvalFunction(MatType);

	if		  (!MatType.compare("GOLD")) { eps = epsAu; mu = muVac;}
	else if (!MatType.compare("SILVER")) { eps = epsAg; mu = muVac; }
	else if (!MatType.compare("ALUMINIUM")) { eps = epsAl; mu = muVac; }
	else if (!MatType.compare("ALUMINUM_NITRIDE")) { eps = epsAlN; mu = muVac; }
	else if (!MatType.compare("CHROMIUM")) { eps = epsCr; mu = muVac; }
	else if (!MatType.compare("TUNGSTEN")) { eps = epsW; mu = muVac; }
	else if (!MatType.compare("PLATINUM")) { eps = epsPt; mu = muVac; }
	else if (!MatType.compare("SILICON_CARBIDE")) { eps = epsSiC; mu = muVac; }
	else if (!MatType.compare("SILICON_DIOXIDE")) { eps = epsSiO2; mu = muVac; }
	else if (!MatType.compare("POTASIUM_BROMIDE")) { eps = epsKBr; mu = muVac; }
	else if (!MatType.compare("ALUMINA")) { eps = epsAl2O3; mu = muVac; }
	else if (!MatType.compare("SILICON_NITRIDE")) { eps = epsSi3N4; mu = muVac; }
	else if (!MatType.compare("IRON_PLATINUM")) { eps = epsFePt; mu = muVac; }
	else if (!MatType.compare("DIAMOND")) { eps = epsC; mu = muVac; }
	else if (!MatType.compare("BORON_NITRIDE")) { eps = epsBN; mu = muVac; }
	else if (!MatType.compare("VACUUM")) { eps = epsVac; mu = muVac; }
	else if (!MatType.compare("PDMS")) { eps = epsPDMS; mu = muVac; }
	else if (!MatType.compare("PMMA")) { eps = epsPMMA; mu = muVac; }
	else if (!MatType.compare("HDPE")) { eps = epsHDPE; mu = muVac; }
	else if (!MatType.compare("LDPE")) { eps = epsLDPE; mu = muVac; }
	else if (!MatType.compare("PVdF-HFP")) { eps = epsPVdF_HFP; mu = muVac; }
	else if (!MatType.compare("VANADIUM_DIOXIDE_R")) { eps = epsVO2h; mu = muVac; }
	else if (!MatType.compare("VANADIUM_DIOXIDE_M")) { eps = epsVO2c; mu = muVac; }
	else if (!MatType.compare("MOLYBDENUM_TRIOXIDE")) { eps = epsMoO3; mu = muVac; }
	else if (!MatType.compare("PERFECT_CONDUCTOR")) { eps0 = DBL_MAX; mu = muVac; }
	else if (!MatType.compare("EPS_FILE")) {
		eps_interp.GetFromFile(Tokens.at(0), false);
		mu = muVac;
	}
	else if (!MatType.compare("CONST_EPS"))
	{
		string epsRe, epsIm;

		// extract real and imaginary part of const_eps
		epsRe = Tokens.at(0);
		epsIm = Tokens.at(1);

		// define value eps0
		eps0 = stod(epsRe) + II*stod(epsIm);
		mu = muVac;
	}

	else { ErrorMsg("Error: No material found"); return 0; }
	return 1;
}

string Material::MaterialType(void) const { return MatType; }

cdouble Material::Eps(double w) {
	// Change units from microns to rad/s (Temporary!!)
	w = 2 * M_PI * SPEEDOFLIGHT / w * 1E6;

	if (!MatType.compare("CONST_EPS")) return eps0;
	if (!MatType.compare("PERFECT_CONDUCTOR")) return eps0;
	if (!MatType.compare("EPS_FILE")) return eps_interp(w);
	return eps(w); }

cdouble Material::Mu(double w) { return mu(w); }

// end of Material Class

// ***************************************************************
// *************** Material Properties Definition
// ***************************************************************

// add logspace option
void printEps(epsmu_fn fx, double w1, double w2, int Nw) {
	double dw = (Nw>1) ? (w2 - w1) / (Nw - 1) : 2*(w2-w1);
	for(double w=w1; w<=w2; w += dw)
		cout << w << " " << fx(w) << endl;
}


cdouble epsAu(double w) {
	//double epsInf = 1., wp = 1.3594122e16, g = 1.0495462e14;
	double epsInf = 1., wp = 1.3713e+16, g = 4.0564e+13;
	return epsInf - wp*wp / (w*w + II*g*w);
}

// for now, just a simple Drude model
// could also have a multi-oscillator model
cdouble epsAg(double w) {
	static cDataInterpol Eps("Ag.txt");
	// create drude model
	double epsInf = 5.0, wp = 1.4889e+16, g = 5.8824e+13;

	if (w > Eps.xmin && w < Eps.xmax) return Eps(w);
	return epsInf - wp * wp / (w*w + II * g*w);
}

// ref: Joulain et. al. PRB 68, 245405 (2003)
cdouble epsAl(double w) {
	double epsInf = 1., wp = 1.747e16, g = 7.596e13;
	return epsInf - wp*wp / (w*w + II*g*w);
}

cdouble epsAlN(double w) {
	static cDataInterpol Eps("AlN.txt");
	return Eps(w);
}

cdouble epsCr(double w) {
	double epsInf = 1., wp = 7.1942e14, g = 7.6921e13;
	return epsInf - wp*wp / (w*w + II*g*w);
}

cdouble epsW(double w) {
	double epsInf = 1., wp = 8.8118e15, g = 7.5963e13;
	return epsInf - wp*wp / (w*w + II*g*w);
}

cdouble epsPt(double w) {
	double wp = 2*M_PI*1.244e15, g = 2*M_PI*16.73e12;
	return 1.0 - wp*wp / (w*w + II*g*w);
}

// SiC, from Ben-Abdallah et. al. J. Appl. Phys. 106, 044306 (2009)
cdouble epsSiC(double w) {
	double epsInf, wLO, wTO, g;
	epsInf = 6.7;
	wLO = 18.253e13;
	wTO = 14.937e13;
	g = 8.966e11;
	return epsInf * (1. + (wLO*wLO - wTO*wTO) / (wTO*wTO - w*w - II*g*w)); // same as in epsCBN
}

cdouble epsSiDoped(double w, double N) {

	double e = ECHARGE;
	double m0 = 9.109e-31;
	double eps0 = EPS0;
	double c=299792458;

	double n =3.4, kappa = 0;
	cdouble epsSi = (n + II*kappa)*(n + II*kappa);
	double m_h = 0.37*m0;         // m_h is effective mass of hole

	// doping concentration and carrier concentration
	double Na=6.9e20*1e6;
	N=N*1e6;

	// parameters for p_type Boron silicon
	double mu_min = 44.9;
	double mu_max = 470.5;
	double Nr = 2.23e17;
	double alpha = 0.719;

	// mu--mobility, unit: cm^2/Vs; mu_min, mu_max, alpha and Nr are fit parameters
	double mu = mu_min+(mu_max-mu_min)/(1+pow(Na/Nr,alpha));

	// convert the unit of 'mu' to m^2/Vs
	mu = mu*1e-4;

	// scattering time
	double tau_h=mu*m_h/e;

    // dielectric function
    return epsSi - N*e*e/(eps0*m_h)/(w*w+II*w/tau_h);
}

// Silica (SiO2)
cdouble epsSiO2(double w) {
	static cDataInterpol Eps("SiO2.txt");
	return Eps(w);
}

// Vanadium Dioxide (ruthile) [VO2 (R)]
cdouble epsVO2h(double w) {
	static cDataInterpol Eps("VO2h.txt");
	return Eps(w);
}

// Vanadium Dioxide (monoclinic) [VO2 (M)]
cdouble epsVO2c(double w) {
	static cDataInterpol Eps("VO2c.txt");
	return Eps(w);
}

// Molybdenum Trioxide (MoO3)
cdouble epsMoO3(double w) {
	static cDataInterpol Eps("MoO3.txt");
	return Eps(w);
}

// Potasium Bromide (KBr)
cdouble epsKBr(double w) {
	static cDataInterpol Eps("KBr.txt");
	return Eps(w);
}

cdouble epsAl2O3(double w) {
	static cDataInterpol Eps("Al2O3.txt");
	return Eps(w);
}

cdouble epsSi3N4(double w) {
	static cDataInterpol Eps("Si3N4.txt");
	return Eps(w);
}

cdouble epsFePt(double w) {
	static cDataInterpol Eps("FePt.txt");
	return Eps(w);
}

// Diamond
cdouble epsC(double w) {
	static cDataInterpol Eps("Diamond.txt");
	return Eps(w);
}

// Polydimethylsiloxane (PDMS)
cdouble epsPDMS(double w) {
	static cDataInterpol Eps("PDMS.txt");
	return Eps(w);
}

// Poly(methyl methacrylate) (PMMA)
cdouble epsPMMA(double w) {
	static cDataInterpol Eps("PDMS.txt");
	return Eps(w);
}

// High density polyethylene (HDPE)
cdouble epsHDPE(double w) {
	static cDataInterpol Eps("HDPE.txt");
	return Eps(w);
}

// Low density polyethylene (LDPE)
cdouble epsLDPE(double w) {
	static cDataInterpol Eps("LDPE.txt");
	return Eps(w);
}

// poly(vinylidene fluoride-co- hexafluoropropene) [P(VdF-HFP)HP].
cdouble epsPVdF_HFP(double w) {
	static cDataInterpol Eps("PVdF-HFP.txt");
	return Eps(w);
}

// c-BN, from Francoeur et. al. JQSRT 110, 2002 (2009)
cdouble epsBN(double w) {
	double epsInf, wLO, wTO, g;
	epsInf = 4.46;
	wLO = 2.451e14;
	wTO = 1.985e14;
	g = 9.934e11;
	return epsInf * (w*w - wLO*wLO + II*g*w) / (w*w - wTO*wTO + II*g*w);
}

cdouble epsVac(double w) {
	return cdouble(1.,0.);
}
cdouble muVac(double w) {
	return cdouble(1., 0.);
}

// for next three methods:
//   with c++11, can use lambda expression to change signature
//   e.g. auto f = [=] (double d) { return epsCont(d,1); };
cdouble epsConst(double w, cdouble eps) {
	static cdouble epslocal = eps;
	return epslocal;
}

cdouble epsDrude(double w, double wp, double g, double epsInf) {
	return cdouble( epsInf - wp*wp /( w*w + g*g ), g*wp*wp /( w*w*w + g*g*w ) );
}
