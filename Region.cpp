#define _USE_MATH_DEFINES
#include "mcphoton_lib.h"

#define EXACT_DIST 0
#define HG_APPROX  1


RTRegion::RTRegion() {
	Initialize();
	Allocate();
}

RTRegion::RTRegion(const string &RegionName,const string &MaterialName)
{
	Initialize();
	Allocate();
	SetRegion(RegionName, MaterialName);
}

RTRegion::~RTRegion()
{
	CleanUp();
}

void RTRegion::Initialize()
{
	Properties = NULL;
	mie_particle = NULL;
	frequency_set = false;
	wExtPath = DBL_MAX;
	is_composite = false;
}

void RTRegion::CleanUp()
{
	if (Properties != NULL) delete Properties;
	if (mie_particle != NULL) {
		for (int i = 0; i < Num_Inclusions; i++){
			mie_particle[i].CleanUp();
		}
		//delete[] mie_particle;
	}
	Initialize();
}

void RTRegion::Allocate()
{
	Properties = new Material;
}

void RTRegion::SetRegion(const string &RegionName,const string &MaterialName)
{
	Label = RegionName;
	OpMaterialName = MaterialName;
	Properties->SetMaterial(MaterialName);
}

string RTRegion::Print(void) const {
	string outstream = "";
	outstream += "	Region Label:	" + Label + "\n";
	outstream += "	Material:	" + Properties->MaterialType() + "\n";
	if (IsComposite())
		for (int i = 0; i < Num_Inclusions; i++)
			outstream += "	Mie Scattering File: " + mie_particle[i].MieFile;
	return outstream;
}

void RTRegion::Set_Inclusions(string strInst) {
	is_composite = true;

	// Tokenize strInst to extract all inclusion declariation
	vector<string> Inst = Tokenize(strInst, "{}"); // tokenize by { }
	for (int i = 0; i < Inst.size(); i++) {
		Inst.at(i) = Inst.at(i).substr(1, Inst.at(i).size() - 2); // remove "{}"
	}
	Num_Inclusions = (int)Inst.size();

	// create values for inclusions and sort index by diameter
	double *np_vol = new double [Num_Inclusions];
	mie_particle = new Mie_Inclusion[Num_Inclusions];
	for (int i = 0; i < Num_Inclusions; i++) {
		mie_particle[i].set_object(Inst.at(i), Label);
		np_vol[i] = (mie_particle[i].NP_Vol);
	}
	sort_xy(np_vol, mie_particle, sizeof(Mie_Inclusion), Num_Inclusions - 1);

	delete[] np_vol;
}

int RTRegion::setProperties(const double &w) {

	if (frequency_set) return 0; // if properties have been already set

	double Deff, lambda, fV;
	cdouble Nh, eps1, mu1;

	// get optical properties of host
	wEps = Properties->Eps(w);
	wMu = Properties->Mu(w);

	// if material is a composite get the scattering properties of inclusions
	if (IsComposite()) {
		for (int i = 0; i < Num_Inclusions; i++){
			if (!mie_particle[i].set_at_frequency(w)) return 0;

			// If particle smaller than wavelength set effective media properties
			if (mie_particle[i].is_opt_properties()) { // if particle has optical properties set
				fV = NP_volfrac_at(i);				// volume fraction of particle
				Deff = NP_EffLength_at(i);			// Effective length of particle (um)
				eps1 = mie_particle[i].get_wEps();	// dielectric constant of particle
				mu1 = mie_particle[i].get_wMu();	// permeability of particle
				Nh = sqrt(wEps*wMu);				// refractive index of host
				lambda = w / Nh.real();				// wavelength in surroundings (um)
				if (Deff < 0.09*lambda) {// if D/lambda << 1 consider particles as effective media
					wEps = eff_Bruggerman(fV, eps1, wEps);
					wMu  = eff_Bruggerman(fV,  mu1,  wMu);
					mie_particle[i].is_scatter = false;		// scattering of particle is not considered anymore
				}
			}
		}
	}

	// get extinction path
	double kext = imag(sqrt(wEps*wMu));
	if (kext > 1E-10) wExtPath = w/ (4 * M_PI * kext); // in um

	frequency_set = true;
	frequency = w;
	return 1;
}

void RTRegion::reset_wProperties(void) {
	wEps = 0;
	wMu = 0;
	if (IsComposite())
		for (int i = 0; i < Num_Inclusions; i++)
			mie_particle[i].reset_wProperties();

	frequency_set = false;
	wExtPath = DBL_MAX;
}

// sample the path length of all the particles and returns the minimum
double RTRegion::NP_sample_path(int &i_np) {
	double *NPlength = new double [Num_Inclusions];
	double *i_particle = new double [Num_Inclusions];

	for (int i = 0; i < Num_Inclusions; i++) {
		if (mie_particle[i].is_scatter)
			NPlength[i] = -log(RandomNum)*mie_particle[i].get_wMeanPath();
		else
			NPlength[i] = DBL_MAX;
		i_particle[i] = i;
	}
	sort_xy(NPlength, i_particle, sizeof(double), Num_Inclusions - 1);

	i_np = (int)i_particle[0]; // return index of particle with shortest path
	double NP_length_min = NPlength[0];
	delete[] NPlength, i_particle;
	return NP_length_min;
}

cdouble RTRegion::get_wEps() {
	return wEps;
}

cdouble RTRegion::get_wMu() {
	return wMu;
}

double RTRegion::get_wExtPath(void) {
	return wExtPath;
}

bool RTRegion::IsComposite(void) const {
	return is_composite;
}

int RTRegion::get_NumInclusions(void) const {
	// Error if region does not contain inclusions
	if (!IsComposite())
		ErrorMsg("Region " + Label + " does not contain inclusions");
	return Num_Inclusions;
}

double RTRegion::NP_concentration_at(const int &i) const {
	// Error if region does not contain inclusions
	if (!IsComposite())
		ErrorMsg("Region " + Label + " does not contain inclusions");

	// Error if the index is larger than the number of inclusions
	if (i >= Num_Inclusions)
		ErrorMsg("Index " + to_string(i) + " larger than the number of inclusions in" + Label);

	// Error if index is negative
	if (i < 0) ErrorMsg("Invalid index at Region " + Label);

	// If pass error return the value
	return mie_particle[i].NPcon;
}

double RTRegion::NP_volfrac_at(const int &i) const {
	// Error if region does not contain inclusions
	if (!IsComposite())
		ErrorMsg("Region " + Label + " does not contain inclusions");

	// Error if the index is larger than the number of inclusions
	if (i >= Num_Inclusions)
		ErrorMsg("Index " + to_string(i) + " larger than the number of inclusions in" + Label);

	// Error if index is negative
	if (i < 0) 
		ErrorMsg("Invalid index at Region " + Label);

	// If pass error return the value
	return mie_particle[i].NPcon*mie_particle[i].NP_Vol;
}

double RTRegion::NP_volume_at(const int &i) const {
	// Error if region does not contain inclusions
	if (!IsComposite())
		ErrorMsg("Region " + Label + " does not contain inclusions");

	// Error if the index is larger than the number of inclusions
	if (i >= Num_Inclusions)
		ErrorMsg("Index " + to_string(i) + " larger than the number of inclusions in" + Label);

	// Error if index is negative
	if (i < 0) 
		ErrorMsg("Invalid index at Region " + Label);

	// If pass error return the value
	return mie_particle[i].NP_Vol;
}

double RTRegion::NP_EffLength_at(const int &i) const {
	// Error if region does not contain inclusions
	if (!IsComposite())
		ErrorMsg("Region " + Label + " does not contain inclusions");

	// Error if the index is larger than the number of inclusions
	if (i >= Num_Inclusions)
		ErrorMsg("Index " + to_string(i) + " larger than the number of inclusions in" + Label);

	// Error if index is negative
	if (i < 0) 
		ErrorMsg("Invalid index at Region " + Label);

	// If pass error return the value
	return mie_particle[i].NP_Leff;
}

string RTRegion::NP_label_at(const int &i) const {
	// Error if region does not contain inclusions
	if (!IsComposite())
		ErrorMsg("Region " + Label + " does not contain inclusions");

	// Error if the index is larger than the number of inclusions
	if (i >= Num_Inclusions)
		ErrorMsg("Index " + to_string(i) + " larger than the number of inclusions in" + Label);

	// Error if index is negative
	if (i < 0) 
		ErrorMsg("Invalid index at Region " + Label);

	// If pass error return the value
	return mie_particle[i].MieFile;
}

int RTRegion::NP_PhaseDistType_at(const int &i) const {
	// Error if region does not contain inclusions
	if (!IsComposite())
		ErrorMsg("Region " + Label + " does not contain inclusions");

	// Error if the index is larger than the number of inclusions
	if (i >= Num_Inclusions)
		ErrorMsg("Index " + to_string(i) + " larger than the number of inclusions in" + Label);

	// Error if index is negative
	if (i < 0) ErrorMsg("Invalid index at Region " + Label);

	// If pass error return the value
	return mie_particle[i].Phase_dist_type;
}

double RTRegion::NP_wMeanPath_at(const int &i) const {
	// Error if region does not contain inclusions
	if (!IsComposite())
		ErrorMsg("Region " + Label + " does not contain inclusions");

	// Error if the index is larger than the number of inclusions
	if (i >= Num_Inclusions)
		ErrorMsg("Index " + to_string(i) + " larger than the number of inclusions in" + Label);

	// Error if index is negative
	if (i < 0) ErrorMsg("Invalid index at Region " + Label);

	// Error if the frequency properties have not been set jet
	if (!frequency_set) ErrorMsg("Properties not set in " + mie_particle[i].MieFile + " in " + Label);

	// If pass error return the value
	return mie_particle[i].get_wMeanPath();
}

double RTRegion::NP_wScaAlbedo_at(const int &i) const {
	// Error if region does not contain inclusions
	if (!IsComposite())
		ErrorMsg("Region " + Label + " does not contain inclusions");

	// Error if the index is larger than the number of inclusions
	if (i >= Num_Inclusions)
		ErrorMsg("Index " + to_string(i) + " larger than the number of inclusions in" + Label);

	// Error if index is negative
	if (i < 0) ErrorMsg("Invalid index at Region " + Label);

	// Error if the frequency properties have not been set jet
	if (!frequency_set) ErrorMsg("Properties not set in " + mie_particle[i].MieFile + " in " + Label);

	// If pass error return the value
	return mie_particle[i].get_wScatAlbedo();
}

double RTRegion::NP_wAsymmetry_at(const int &i) const {
	// Error if region does not contain inclusions
	if (!IsComposite())
		ErrorMsg("Region " + Label + " does not contain inclusions");

	// Error if the index is larger than the number of inclusions
	if (i >= Num_Inclusions)
		ErrorMsg("Index " + to_string(i) + " larger than the number of inclusions in" + Label);

	// Error if index is negative
	if (i < 0) ErrorMsg("Invalid index at Region " + Label);

	// Error if the frequency properties have not been set jet
	if (!frequency_set) ErrorMsg("Properties not set in " + mie_particle[i].MieFile + " in " + Label);

	// If pass error return the value
	return mie_particle[i].get_wAsym();
}

dataInterpol RTRegion::NP_wfDist_at(const int &i) const {
	// Error if region does not contain inclusions
	if (!IsComposite())
		ErrorMsg("Region " + Label + " does not contain inclusions");

	// Error if the index is larger than the number of inclusions
	if (i >= Num_Inclusions)
		ErrorMsg("Index " + to_string(i) + " larger than the number of inclusions in" + Label);

	// Error if index is negative
	if (i < 0) ErrorMsg("Invalid index at Region " + Label);

	// Error if the frequency properties have not been set jet
	if (!frequency_set) ErrorMsg("Properties not set in " + mie_particle[i].MieFile + " in " + Label);

	// If pass error return the value
	return mie_particle[i].get_wfDist();
}

cdouble RTRegion::NP_wEps_at(const int &i) const {
	// Error if region does not contain inclusions
	if (!IsComposite())
		ErrorMsg("Region " + Label + " does not contain inclusions");

	// Error if the index is larger than the number of inclusions
	if (i >= Num_Inclusions)
		ErrorMsg("Index " + to_string(i) + " larger than the number of inclusions in" + Label);

	// Error if index is negative
	if (i < 0) ErrorMsg("Invalid index at Region " + Label);

	// Error if the frequency properties have not been set jet
	if (!frequency_set) ErrorMsg("Properties not set in " + mie_particle[i].MieFile + " in " + Label);

	// If pass error return the value
	return mie_particle[i].get_wEps();
}

cdouble RTRegion::NP_wMu_at(const int &i) const {
	// Error if region does not contain inclusions
	if (!IsComposite())
		ErrorMsg("Region " + Label + " does not contain inclusions");

	// Error if the index is larger than the number of inclusions
	if (i >= Num_Inclusions)
		ErrorMsg("Index " + to_string(i) + " larger than the number of inclusions in" + Label);

	// Error if index is negative
	if (i < 0) ErrorMsg("Invalid index at Region " + Label);

	// Error if the frequency properties have not been set jet
	if (!frequency_set) ErrorMsg("Properties not set in " + mie_particle[i].MieFile + " in " + Label);

	// If pass error return the value
	return mie_particle[i].get_wMu();
}
