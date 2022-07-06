#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <math.h>
#include "mcphoton_lib.h"
//#include "localmathlib.h"

Esource::Esource(int Nhw, 
	string wfunc, string phifun, string thetafun, string Edistfun) {
	Initiallize();
	Nphotons = Nhw;
	setWavelength(wfunc);
	setOrientation(phifun, thetafun);
	setEfield(Edistfun);
}

Esource::~Esource(void) {
	CleanUp();
}

void Esource::Initiallize(void) {
	w = NULL;
	theta = NULL;
	phi = NULL;
	Nphotons = 0;
	polarization = "RANDOM";
	Epower.Val = 0;
	Epower.strfunc = "POINT_SOURCE";
	iPhi = 0, iThe = 0, iw = -1;
	next_step = true;
}

void Esource::CleanUp(void) {
	if (w     != NULL) delete[] w;
	if (theta != NULL) delete[] theta;
	if (phi   != NULL) delete[] phi;
	Initiallize();
}

void Esource::setWavelength(string func) {
	vstring token;
	vstring arg;

	token = Tokenize(func);
	wUnit = token.at(1);

	func = token.at(0);
	w = SetRange(func,Nw);
	Convert_wUnits();
}

void Esource::setOrientation(string strPhi, string strTheta) {
	theta = SetAngle(strTheta, Ntheta, thetaUnit);
	phi = SetAngle(strPhi, Nphi, phiUnit);
}

void Esource::setEfield(string strInst) {
	string EfieldID = "Photon Source";
	objArgsearch EfielDist[] = {
		{ "CENTER"		, PA_POINT3D , false, 1, 1, (void *)&Beam_center  },
		{ "POLARIZATION", PA_STRING  , false, 1, 1, (void *)&polarization },
		{ "POWER"       , PA_DBLSTR  , false, 1, 1, (void *)&Epower       },
		{ "",0,0,0,0,0 }
	};

	GetArguments(strInst, EfielDist, 3, EfieldID);
}

Photon* Esource::MakePhoton(RTRegion **Regions, int NumRegions) const{
	double w_hw = w[iw];
	Photon *hw = new Photon(w_hw);
	Point3D  r0 = GeneratePosition(theta[iThe], phi[iPhi]);
	Vector3D hk = GenerateMomentum(theta[iThe], phi[iPhi]);
	Vector3D Epol = GeneratePolarization();

	// pass data to photon
	hw->SetPosition(r0);
	hw->SetMomentum(hk);
	hw->SetPolarization(Epol);

	// Set the region of the photon
	for (int ir = 0; ir < NumRegions; ir++) {
		if (!Regions[ir]->Label.compare("EXTERIOR")) {
			hw->SetRegion(Regions[ir]);
		}
	}
	return hw;
}

Point3D Esource::GeneratePosition(double xtheta, double xphi) const{
	double r, Tol = 1E-10;
	Point3D pos;

	string Edist = Epower.strfunc;
	double rDist = Epower.Val;
	// get radial possition of photon according to distribution
	if        (!Edist.compare("POINT_SOURCE")) {
		return Beam_center;
	}
	
	if (!Edist.compare("RADIAL_UNIFORM")) {
		r = RandomNum*rDist;
	} else if (!Edist.compare("GAUSSIAN")) {
		r = rDist*sqrt(-log(1 - RandomNum + 1E-10));
	}

	// generate possition in the xy plane
	if (r > 0){
		double x1, y1;
		for (;;) {	/*new direction*/
			x1 = 2.0*RandomNum - 1.0; //sample cos(phi)
			y1 = 2.0*RandomNum - 1.0; //sample sin(phi)
			if (x1*x1 + y1*y1 == 1) break;
		}
		pos.SetPoint3D(x1*r, y1*r, 0); // possition in the plane of the beam
	} else {
		pos.SetPoint3D(0, 0, 0); // possition in the plane of the beam
	}

	// pass local position to global coordinates
	xtheta = M_PI - xtheta;
	double cosT = cos(xtheta);
	double sinT = sqrt(1 - cosT*cosT);
	double cosP = sinT < Tol ? cosT / abs(cosT) : cos(xphi);
	double sinP = sinT < Tol ? 0 : sin(xphi);
	
	Point3D rxx;
	rxx.x = +cosP*cosT*pos.x - sinP*pos.y + cosP*sinT*pos.z;
	rxx.y = +sinP*cosT*pos.x + cosP*pos.y + sinP*sinT*pos.z;
	rxx.z = -sinT*pos.x + 0 * pos.y + cosT*pos.z;

	return rxx + Beam_center; // return photon's possition
}

Vector3D Esource::GeneratePolarization(void) const{
	double pol_x, pol_y;

	if (!polarization.compare("RANDOM")) {
		double psi = 2.0*M_PI*RandomNum;
		pol_x = cos(psi);
		pol_y = sin(psi);	/* sintheta is always positive */
	}
	else if(!polarization.compare("S")) {
		pol_x = 0;
		pol_y = 1;
	}
	else if (!polarization.compare("P")) {
		pol_x = 1;
		pol_y = 0;
	}
	return Vector3D(pol_x, pol_y);
}

Vector3D Esource::GenerateMomentum(double theta, double phi) const{
	double hk_x, hk_y, hk_z;
	hk_x = sin(M_PI - theta)*cos(phi);
	hk_y = sin(M_PI - theta)*sin(phi);
	hk_z = cos(M_PI - theta);
	return Vector3D(hk_x, hk_y, hk_z);
}

bool Esource::RunSimulation(void) {
	iw++;
	if (iw == Nw) { iw = 0; iThe++; }
	if (iThe == Ntheta) { iThe = 0; iPhi++; }
	if (iPhi == Nphi) { next_step = false; }
	return next_step;
}

bool Esource::End_spectrum(void) const{

	if (iw == Nw - 1) return true;
	return false;
}

bool Esource::Start_spectrum(void) const {

	if (iw == 0) return true;
	return false;
}

double *Esource::SetAngle(string str, int &Nx, string &Unit) {
	vstring token;
	vstring arg;
	double xC;

	token = Tokenize(str);
	Unit = token.at(1);
	if       (!Unit.compare("deg")) xC = M_PI / 180.0;
	else if  (!Unit.compare("rad")) xC = 1;
	else {
		ErrorMsg("Esource: Angle unit not defined");
	}
	return SetRange(token.at(0),Nx, xC);
}

double *Esource::SetRange(string func, int &Nx, double C) {
	vstring token;
	vstring arg;
	double Xmin, Xmax;

	arg = EvalFunction(func);
	if (arg.empty()) {
		Xmin = (double)stod(func); Xmin *= C;
		Xmax = (double)stod(func); Xmax *= C;
		Nx = 1;
		return linspace(Xmin, Xmax, Nx);
	}
	else if (!func.compare("linspace")) {
		Xmin = (double)stod(arg.at(0)); Xmin *= C;
		Xmax = (double)stod(arg.at(1)); Xmax *= C;
		Nx   = (int)stod(arg.at(2));
		return linspace(Xmin, Xmax, Nx);
	}
	else if (!func.compare("logspace")) {
		Xmin = (double)stod(arg.at(0)); Xmin *= C;
		Xmax = (double)stod(arg.at(1)); Xmax *= C;
		Nx   = (int)stod(arg.at(2));
		return logspace(Xmin, Xmax, Nx);
	}

	return 0;
}

/*****************************************************************************************/
// Frequency data out
int Esource::GetWavelength_idx(void) const { return iw; }
int Esource::GetWavelength_N(void) const { return Nw; }

double Esource::GetWavelength_val(bool printout) const {
	 // if printout false retrieve frequency in microns
	if (!printout) return w[iw];

	double w_out;
	if      (!wUnit.compare("rad/s")) w_out = 2 * M_PI*SPEEDOFLIGHT / w[iw] * 1E6;
	else if (!wUnit.compare("Hz")) w_out = SPEEDOFLIGHT / w[iw] * 1E6;
	else if (!wUnit.compare("eV")) w_out = 2 * M_PI*SPEEDOFLIGHT / (EV*w[iw]/HBAR) * 1E6;
	else if (!wUnit.compare("nm")) w_out = w[iw] * 1E-3;
	else if (!wUnit.compare("um")) w_out = w[iw];
	else if (!wUnit.compare("mm")) w_out = w[iw] * 1E3;
	else ErrorMsg("Esource: Frequency unit not recognized");

		return w_out;
}

double* Esource::GetWavelength_array(bool printout) const {
	if (!printout) return w;
	
	double *w_out = new double[Nw];
	for (int iiw = 0; iiw < Nw; iiw++) {
		if      (!wUnit.compare("rad/s")) w[iiw] = 2 * M_PI*SPEEDOFLIGHT / w[iiw] * 1E6;
		else if (!wUnit.compare("Hz"))  w[iiw] = SPEEDOFLIGHT / w[iiw] * 1E6;
		else if (!wUnit.compare("eV")) w[iiw] = 2 * M_PI*SPEEDOFLIGHT / (EV*w[iiw]/HBAR)*1E6;
		else if (!wUnit.compare("nm")) w[iiw] = w[iiw] * 1E-3;
		else if (!wUnit.compare("um")) {}
		else if (!wUnit.compare("mm")) w[iiw] = w[iiw] * 1E+3;
		else  ErrorMsg("Esource: Frequency unit not recognized");
	}
	return w_out;
}
/*****************************************************************************************/


/*****************************************************************************************/
// Zenith incidence angle data out
double* Esource::GetZenith_array(void) const { return theta;}
int Esource::GetZenith_N(void) const { return Ntheta; }
int Esource::GetZenith_idx(void) const { return iThe; }
double Esource::GetZenith_val(bool printout) const {
	if (!printout) { // if printout false retrieve frequency in rad/s
		return theta[iThe];;
	}
	else {
		double theta_out;
		if (!thetaUnit.compare("rad")) { theta_out = theta[iThe]; }
		else if (!thetaUnit.compare("deg")) { theta_out = theta[iThe] * 180.0/ M_PI; }
		else { ErrorMsg("Esource:theta unit not recognized"); }

		return theta_out;
	}
}
/*****************************************************************************************/

/*****************************************************************************************/
// Azimuth incidence angle data out
double* Esource::GetAzimuth_array(void) const { return phi; }
int Esource::GetAzimuth_N(void) const { return Nphi; }
int Esource::GetAzimuth_idx(void) const { return iPhi; }
double Esource::GetAzimuth_val(bool printout) const {
	if (!printout) { // if printout false retrieve frequency in rad/s
		return phi[iPhi];;
	}
	else {
		double phi_out;
		if (!phiUnit.compare("rad")) { phi_out = phi[iPhi]; }
		else if (!phiUnit.compare("deg")) { phi_out = phi[iPhi] * 180.0 / M_PI; }
		else { ErrorMsg("Esource: phi unit not recognized"); }

		return phi_out;
	}
}
/*****************************************************************************************/

int Esource::GetTotalPhotons() const {
	return Nphotons;
}

bool Esource::FirstRun() const {
	if (iw == 0 && iThe == 0 && iPhi == 0) return true;
	else return false;
}

// convert the input wavelength into microns
void Esource::Convert_wUnits() {
	for (int iiw = 0; iiw < Nw; iiw++){
		if      (!wUnit.compare("rad/s")) w[iiw] = 2 * M_PI*SPEEDOFLIGHT / w[iiw] * 1E6;
		else if (!wUnit.compare("Hz")) w[iiw] = SPEEDOFLIGHT / w[iiw] * 1E6;
		else if (!wUnit.compare("eV")) w[iiw] = SPEEDOFLIGHT / (EV * w[iiw] / HBAR) * 1E6;
		else if (!wUnit.compare("nm")) w[iiw] = w[iiw] * 1E-3;
		else if (!wUnit.compare("um")) {}
		else if (!wUnit.compare("mm")) w[iiw] = w[iiw] * 1E+3;
		else  ErrorMsg("Esource: Frequency unit not recognized");
	}

}