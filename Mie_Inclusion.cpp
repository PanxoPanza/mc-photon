#define _USE_MATH_DEFINES
#include "mcphoton_lib.h"

#define EXACT_DIST 0
#define HG_APPROX  1

Mie_Inclusion::Mie_Inclusion(void) {
	Initialize();
}

Mie_Inclusion::Mie_Inclusion(const string &str_inst, const string &xRegionLabel) {
	Initialize();
	set_object(str_inst, xRegionLabel);
}

Mie_Inclusion::~Mie_Inclusion(void) {
	CleanUp();
}

void Mie_Inclusion::set_object(const string &str_inst, const string &xRegionLabel) {
	string MieScattFile[2];
	MieScattFile[1] = "HG_APPROX"; // Henyey-Greenstein approximation taken by default
	vector<string> args;
	bool checkVolFrac = false, checkMieFile = false;

	// set region label
	RegionLabel = xRegionLabel;
	string RegionID = "Region " + RegionLabel;

	// get arguments from the string input "str_inst"
	objArgsearch Inclusion[] = {
		{"VOLFRAC", PA_DOUBLE, false, 1, 1, (void *)&fV},
		{"MIEFILE", PA_STRING, false, 1, 2, (void *)&MieScattFile},
		{"",0,0,0,0,0}
	};
	GetArguments(str_inst, Inclusion, 2, RegionID);

	// Store filename
	MieFile = MieScattFile[0];

	// Get distribution function type
	if (!MieScattFile[1].compare("EXACT"))			Phase_dist_type = EXACT_DIST;
	else if (!MieScattFile[1].compare("HG_APPROX"))	Phase_dist_type = HG_APPROX;
	else ErrorMsg("Undefined phase distribution in: " + MieFile + " from " + RegionLabel);

	// if the phase distribution data is not available set HG_APPROX
	if (!is_pDist_exat) Phase_dist_type = HG_APPROX; // PENDING! send a message to user

	// get data from *.mie file
	Open_MieFile(MieScattFile);

}

void Mie_Inclusion::Initialize(void) {

	fV = 0;
	NP_Leff = 0;
	NP_Vol = 0;
	NPcon = 0;
	NP_Label.clear();
	Host_Label.clear();
	Mie_calculate = false;
	NumFreq = 0;

	MieFileOpen = false;
	MieFile = "";
	RegionLabel = "";
	cosT.clear();

	MieData = false;
	sCross.clear();
	aCross.clear();
	gAsymm.clear();

	eff_opt = false;
	eps[0].clear(); eps[1].clear();
	mu[0].clear();  mu[1].clear();

	is_pDist_exat = false;
	Ntheta = 0;
	cosT.clear();
	fTheta = NULL;

	reset_wProperties();
}

void Mie_Inclusion::reset_wProperties(void) {
	frequency_set = false;
	is_scatter = true;
	wMeanPathParticle = 0;
	wScatAlbedo = 0;
	wgAsym = 0;
	wEps = 0;
	wMu = 0;
	wNPftheta.clear();
}

void Mie_Inclusion::CleanUp(void) {
	Initialize();
	if (fTheta != NULL) delete[] fTheta;
}

void Mie_Inclusion::Open_MieFile(string MieFileName[]) {

	ifstream miescatData(MieFileName[0] + ".mie");
	if (!miescatData.is_open())
		ErrorMsg("Unable to open: " + MieFileName[0] + ".mie");

	// File successfuly opened
	MieFileOpen = true;
	if (!Read_Parameters(miescatData))
		ErrorMsg("Not enough parameters in " + MieFile + ".mie");
	if (!Read_Data(miescatData))
		ErrorMsg("Failed to read data in " + MieFile + ".mie");
	miescatData.close();

	// Get particle's concentration
	NPcon = fV / NP_Vol;
}

int Mie_Inclusion::Read_Parameters(ifstream &File) {
	string line;
	vector<string> token;
	string begobj = "MIEPARAM" , endobj = "END_MIEPARAM";
	bool objstate = false;

	objArgsearch Inclusion[] = {
	{"NP_LABEL"			, PA_STRING	, false, 1, 1, (void *)&NP_Label	},
	{"HOST_LABEL"		, PA_STRING	, false, 1, 1, (void *)&Host_Label	},
	{"NP_VOLUME"	    , PA_DOUBLE	, false, 1, 1, (void *)&NP_Vol	    },
	{"NP_LENGTH"	    , PA_DOUBLE	, false, 1, 1, (void *)&NP_Leff	    },
	{"NUM_FREQ"			, PA_INT	, false, 1, 1, (void *)&NumFreq		},
	{"PHASEDIST_ANGLE"	, PA_INT	, false, 1, 1, (void *)&Ntheta		},
	{"MIE_DATA"			, PA_BOOL	, false, 1, 1, (void *)&MieData		},
	{"OPT_PROPERTIES"	, PA_BOOL	, false, 1, 1, (void *)&eff_opt		},
	{"",0,0,0,0,0}
	};
	while (getline(File, line)) {
		CommentOut(line); // remove comments
		token = Tokenize(line); // tokenize considering " " - "\t" - "," characters

		// End of object declaration
		if (!token.empty()) { // if line is not empty
			if (objstate && !token.front().compare(endobj)) {
				objstate = false; // End of object declaration
				break;
			}

			// Object found (extract arguments)
			if (objstate) GetArguments(token.front(), Inclusion, 1, MieFile);
			if (!token.front().compare(begobj)) objstate = true; // Start of object declaration
		}
	}

	File.clear();
	File.seekg(0, ios::beg);
	if (Ntheta > 0) is_pDist_exat = true;

	// check all input paramters
	if (NP_Label.empty()) return 0;					// if no particle label
	if (Host_Label.empty()) return 0;				// if no host label
	if (NP_Vol <= 0) return 0;					    // if volume is <= 0
	if (MieData && NumFreq <= 0) return 0;			// if mie-data is true but not frequency range
	if (eff_opt && NumFreq <= 0) return 0;			// if eff. opt. properties is true but not frequency range
	if (is_pDist_exat && NumFreq <= 0) return 0;	// if phase dist. is true but not frequency range

	// **Future project: This is by giving the code the ability to compute mie-scattering
	// if no mie-data and no material properties have been given return 0
	if (!MieData && !Mie_calculate) return 0;

	return 1;
}

int Mie_Inclusion::Read_Data(ifstream &File) {
	string line;
	vector<string> token;
	string begobj = "MIEDATA", endobj = "END_MIEDATA";
	bool objstate = false;
	int id;

	// Read data from file
	dvector wData, sData, aData, gData; // MieData
	dvector Eps_data[2], Mu_data[2]; // Effective properties
	dvector *ftData;

	// Extract scattering phase function
	if (is_pDist_exat) ftData = new dvector[Ntheta];

	while (getline(File, line)) {
		CommentOut(line); // remove comments
		token = Tokenize(line); // tokenize considering " " - "\t" - "," characters

		// End of object declaration
		if (!token.empty()) { // if line is not empty
			if (objstate && !token.front().compare(endobj)) {
				objstate = false; // End of object declaration
				break;
			}

			// Object found (extract arguments)
			if (objstate) {
				id = 0;
				if (MieData) {
					if (id < token.size()) wData.push_back(stod(token.at(id++)));
					if (id < token.size()) sData.push_back(stod(token.at(id++)));
					if (id < token.size()) aData.push_back(stod(token.at(id++)));
					if (id < token.size()) gData.push_back(stod(token.at(id++)));
				}
				if (eff_opt) {
					if (id < token.size()) Eps_data[0].push_back(stod(token.at(id++)));
					if (id < token.size()) Eps_data[1].push_back(stod(token.at(id++)));
					if (id < token.size()) Mu_data[0].push_back(stod(token.at(id++)));
					if (id < token.size()) Mu_data[1].push_back(stod(token.at(id++)));
				}
				if (is_pDist_exat) {
					for (int i = 0; i < Ntheta; i++)
						if (id < token.size()) ftData[i].push_back(stod(token.at(id++)));
				}
			}
			if (!token.front().compare(begobj)) objstate = true; // Start of object declaration
		}
	}
	// check if extracted data does not have any errors
	for (int i = 0; i < sData.size(); i++)
		if (sData.at(i) < 0)
			ErrorMsg("Error in" + MieFile + ".mie: " + "Negative scat. data at row: " + to_string(i + 1));

	for (int i = 0; i < aData.size(); i++)
		if (aData.at(i) < 0)
			ErrorMsg("Error in" + MieFile + ".mie: " + "Negative abs. data at row: " + to_string(i + 1));

	// Correct scattering considering interparticle distance (pending)
	/*
	double cinter = (0.905 / pow(fV, 1.0 / 3.0) - 1)*NP_Leff; // Interparticle's separation
	for (int i = 0; i < sData.size(); i++){
		double lambda = 2 * M_PI*SPEEDOFLIGHT / wData.at(i)*1E6;
		sData.at(i) = sData.at(i)*pow(10, -pow(10, 0.25 - 5.1*cinter / lambda));
	}
	*/

	// Pass data to interpolation function
	if (MieData) {
		sCross.set_points(wData, sData); sData.clear();
		aCross.set_points(wData, aData); aData.clear();
		gAsymm.set_points(wData, gData); gData.clear();
	}
	if (eff_opt) {
		eps[0].set_points(wData, Eps_data[0]); Eps_data[0].clear();
		eps[1].set_points(wData, Eps_data[1]); Eps_data[1].clear();
		mu[0].set_points(wData, Mu_data[0]); Mu_data[0].clear();
		mu[1].set_points(wData, Mu_data[1]); Mu_data[1].clear();
	}
	if (is_pDist_exat) {
		fTheta = new dataInterpol [Ntheta];
		double *xcosT = linspace(-1, 1, Ntheta);
		for (int it = 0; it < Ntheta; it++) {
			cosT.push_back(xcosT[it]);
			fTheta[it].set_points(wData, ftData[it]); 
			ftData[it].clear();
		}
		delete[] ftData;
		delete[] xcosT;
	}
	wData.clear();
	return 1;
}

int Mie_Inclusion::setNPfTheta(const double &w) {
	dvector vfTheta;
	if (frequency_set) return 0; // if properties have been already set

	for (int it = 0; it < Ntheta; it++) {
		vfTheta.push_back(fTheta[it](w));
	}
	wNPftheta.set_points(vfTheta, cosT);
	return 1;
}

int Mie_Inclusion::set_at_frequency(const double &w) {

	if (frequency_set) return 0; // if properties have been already set

	// Mean path at w
	double wDecay = NPcon*(aCross(w) + sCross(w));
	wMeanPathParticle = wDecay < 1E-20? DBL_MAX: 1.0 / wDecay;

	// Scattering albedo  = Csca / (Csca + Cabs), at w
	wScatAlbedo = sCross(w) / (aCross(w) + sCross(w));

	// Optical properties at w
	if (eff_opt){
		wEps = eps[0](w) + II*eps[1](w);
		wMu  =  mu[0](w) +  II*mu[1](w);
	}

	switch (Phase_dist_type) {
	case EXACT_DIST:
		if (!setNPfTheta(w)) return 0;
		break;
	case HG_APPROX:
		wgAsym = gAsymm(w);
		break;
	default:
		ErrorMsg("Unknown phase distribution function in: " + MieFile + " from " + RegionLabel);
	}
	frequency_set = true;

	/*
	Log("NPcon: %.4f",NPcon);
	Log("mu_a (mm^-1): %.4f",NPcon*aCross(w));
	Log("mu_s (mm^-1): %.4f",NPcon*sCross(w));
	Log("wMeanPathParticle: %.4f",wMeanPathParticle);
	Log("wScatAlbedo: %.4f",wScatAlbedo);
	Log("wgAsym: %.4f",wgAsym);
	*/
	
	return 1;
}

bool Mie_Inclusion::wProperties_set(void) const {
	return frequency_set;
}

double Mie_Inclusion::get_wMeanPath(void)const {
	if (!frequency_set)
		ErrorMsg("Properties not set for " + NP_Label + " in " + RegionLabel);

	return wMeanPathParticle;
}

double Mie_Inclusion::get_wScatAlbedo(void)const {
	if (!frequency_set)
		ErrorMsg("Properties not set for " + NP_Label + " in " + RegionLabel);

	return wScatAlbedo;
}

double Mie_Inclusion::get_wAsym(void)const {
	if (!frequency_set)
		ErrorMsg("Properties not set for " + NP_Label + " in " + RegionLabel);

	return wgAsym;
}

cdouble Mie_Inclusion::get_wEps(void)const {
	if (!frequency_set)
		ErrorMsg("Properties not set for " + NP_Label + " in " + RegionLabel);

	return wEps;
}

cdouble Mie_Inclusion::get_wMu(void)const {
	if (!frequency_set)
		ErrorMsg("Properties not set for " + NP_Label + " in " + RegionLabel);

	return wMu;
}

dataInterpol Mie_Inclusion::get_wfDist(void)const {
	if (!frequency_set)
		ErrorMsg("Properties not set for " + NP_Label + " in " + RegionLabel);

	return wNPftheta;
}

bool Mie_Inclusion::is_opt_properties(void) const {
	return eff_opt;
}

bool Mie_Inclusion::is_mie_data(void) const {
	return MieData;
}

bool Mie_Inclusion::is_phasedist(void) const {
	return is_pDist_exat;
}
