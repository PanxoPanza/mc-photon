#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include "MCRT_library.h"


// plot types for monitor
#define RAD_PROP 0
#define FRESNEL  1

/**************************************** Subclass RTSurface ****************************************/
RTSurface::RTSurface() : Surface() {}

RTSurface::RTSurface(RTRegion *External, RTRegion *Internal, string &GeometryArg, const string &xLabel, const int &idx) :
	Surface(External, Internal, GeometryArg, xLabel, idx) {
	NumLayers = 0;
	NumRegions = 2;
	matFilm = NULL;
	tFilm = NULL;
	epsFilm = NULL;
	muFilm = NULL;
	frequency_set = false;
}

RTSurface::~RTSurface(){ 
	CleanUp();

	if (matFilm != NULL) delete[] matFilm;
	if (tFilm   != NULL) delete[] tFilm;
	if (epsFilm != NULL) delete[] epsFilm;
	if (muFilm  != NULL) delete[] muFilm;
	NumLayers = 0;
	matFilm = NULL;
	tFilm = NULL;
	epsFilm = NULL;
	muFilm = NULL;
}

void RTSurface::SetCoating(string strThinFilms) {
	doublestr LayerInfo;
	vector<string> token;
	string SurfaceID = "Surface (" + to_string(Index) + ")";
	objArgsearch SurfSet[] = {
		{"LAYER"			, PA_DBLSTR , false, 1, 1, (void *)&LayerInfo	},
		{"",0,0,0,0}
	};

	token = Tokenize(strThinFilms);// tokenize layer functions
	NumLayers = (int)token.size();
	matFilm = new Material[NumLayers];
	tFilm = new double[NumLayers];
	epsFilm = new cdouble[NumLayers];
	muFilm = new cdouble[NumLayers];

	for (int i = 0; i < NumLayers; i++) {
		GetArguments(token.at(i), SurfSet, 1, SurfaceID);// Extract parameters of each layer
		tFilm[i] = LayerInfo.Val*1E-6;
		matFilm[i].SetMaterial(LayerInfo.strfunc);
		SurfSet->IsFound = false;// Set to false to get next layer
	}
}

bool RTSurface::PhotonHits(Photon &hw) const {
	int iPanel;
	bool IsHit = false;
	Point3D x0, x0_new;
	double dz = DBL_MAX, dz_new;

	for (int i = 0; i < NumPanels; i++) {
		if (Panels[i].PhotonHits(hw, dz_new, x0_new)){
			if (dz_new < dz) {
				IsHit = true;
				iPanel = i;
				dz = dz_new;
				x0 = x0_new;
			}
		}
	}
	if (IsHit) hw.SetSurfaceHit(x0, dz, iPanel, Index);
	return IsHit;
}

int RTSurface::RefractPhoton(Photon &hw) const {
	int iFrom, iTo, kz_Dir;
	if (hw.PhotonIs(Regions[0])) {
		iFrom = 0;  // index of region where photon comes
		iTo = 1;    // index of region where photon goes
		kz_Dir = 1;
	}
	else if (hw.PhotonIs(Regions[1])) {
		iFrom = 1;  // index of region where photon comes
		iTo = 0;    // index of region where photon goes
		kz_Dir = -1;
	}
	else {
		ErrorMsg("Error in Surface: " + Label + ". Photon is neither in Region 1 or 2");
		return 2;
	}

	int FresnelEvent = FresnelRefraction(iFrom, iTo, kz_Dir, hw);

	if (FresnelEvent == 1) { // Photon reflected. Returns to "From" region
		hw.SetRegion(Regions[iFrom]);
	}
	else if (FresnelEvent == -1) { // Photon transmitted to a new region
		hw.SetRegion(Regions[iTo]);
	}
	else if (FresnelEvent == 0) { // Photon absorbed at the interface - Kill
		hw.PhotonKill();
	}
	//cout << hw.Print();
	return FresnelEvent;
}

int RTSurface::set_wLayers(const double &w) {

	if (frequency_set) return 0; // if properties have been already set

	if (NumLayers > 0) {
		for (int i = 0; i < NumLayers; i++) {
			epsFilm[i] = matFilm[i].Eps(w);
			muFilm[i] = matFilm[i].Mu(w);
		}
	}

	frequency_set = true;
	return 1;
}

void RTSurface::reset_wProperties(void) {

	if (NumLayers > 0) {
		for (int i = 0; i < NumLayers; i++) {
			epsFilm[i] = cdouble(0.0, 0.0);
			muFilm[i] = cdouble(0.0, 0.0);
		}
	}
	frequency_set = false;
}


void RTSurface::TransferMatrix(const cdouble &kx0, double w, int kz_Dir, string pol,
	cdouble &m11, cdouble &m12, cdouble &m21, cdouble &m22) const {

	int i;
	double c0 = SPEEDOFLIGHT;// d;
	cdouble eps, mu, Ni, pf, kzd, ki;
	cdouble m11_new, m12_new, m21_new, m22_new;
	cdouble sinTi, cosTi;
	double d;
	for (int j = 0; j < NumLayers; j++) {
		switch (kz_Dir) {
		case  1:  
			i = j;			// wave propagate in the possitive z-direction
			break;
		case -1:  
			i = NumLayers - 1 - j;// wave propagate in the negative  z-direction
			break;
		}
		
		// Set properties of the film
		eps = epsFilm[i]; // dielectric constant of the film
		mu  = muFilm[i];   // relative permeability of the film
		d = tFilm[i]; // film thickness
		Ni = sqrt(eps*mu) != 0.0? sqrt(eps*mu) : 1E-15;
		ki = w*Ni / c0; // wavevector modulus in the film
		sinTi = kx0 / ki;
		cosTi = sqrt(1.0 - sinTi*sinTi); // kz component in the film

		cosTi = cosTi != 0.0 ? cosTi : 1E-15; // make cosTi a small number if zero
		sinTi = sinTi != 0.0 ? sinTi : 1E-15; // make sinTi a small number if zero

		kzd = ki * cosTi* d;
		if (pol.compare("s") == 0) {
			pf = - Z0*mu/(Ni*cosTi);// surface impedance
		}
		else if (pol.compare("p") == 0) {
			pf =   Z0*mu*cosTi/Ni; // surface impedance
		}

		// Get transfer matrix coeffients
		m11_new =        cos(kzd)*m11 - II/pf*sin(kzd)*m12;
		m12_new = -II*pf*sin(kzd)*m11 +       cos(kzd)*m12;
		m21_new =        cos(kzd)*m21 - II/pf*sin(kzd)*m22;
		m22_new = -II*pf*sin(kzd)*m21 +       cos(kzd)*m22;

		// store new values
		m11 = m11_new;
		m12 = m12_new;
		m21 = m21_new;
		m22 = m22_new;
	}
}

int RTSurface::FresnelRefraction(const int &idx_from, const int &idx_to, const int &kz_Dir, Photon &hw) const
{
	/*************************************************************************************
	Refraction of a photon at a planar interface between medium 1 and 2
	based on Fresnel's Law for a EM-wave coming from medium 1.

	Francisco Ramirez 10/2018
	--------------------------------------------------------------------------------
	This function performs refraction from Fresnel coefficients. The optical properties
	of the origin (eps1 and mu1) are obtained from the photon's information (hw.GetRegion),
	while the destiny's optical properties (eps2 and mu2) are set from the calling Surface
	object.
	**************************************************************************************/

	int FresnelState; // 1 if photon reflected, -1 if transmitted, 0 if absorbed
	int PanelID;
	double Tol = 1E-10;
	cdouble p1, p2, ct;
	cdouble m11, m12, m21, m22;
	Point3D x0_loc = hw.x_SurfaceHits();

	// Select the panel hited by the photon
	if (hw.GetSurfacePanelID() >= 0) {
		PanelID = hw.GetSurfacePanelID();
	}
	else {
		ErrorMsg("Error in Surface(" + to_string(Index) + "): PanelID has not being set");
		return 3;
	}

	// Get input parameters
	double w = hw.GetFrequency(); // rad/s
   
    //optical properties medium 1
	cdouble eps1 = Regions[idx_from]->get_wEps(); 
	cdouble mu1  = Regions[idx_from]->get_wMu();
	cdouble n1   = sqrt(eps1*mu1) =! 0.0 ? sqrt(eps1*mu1) : 1E-15;

	//optical properties medium 2
	cdouble eps2 = Regions[idx_to]->get_wEps();
	cdouble mu2  = Regions[idx_to]->get_wMu();
	cdouble n2   = sqrt(eps2*mu2) =! 0.0 ? sqrt(eps2*mu2) : 1E-15;

	// Get vectors in panel's local coordinates
	Vector3D k_hat = Panels[PanelID].GlobaltoLocal(hw.k_hat()); // photon momentum
	Vector3D e_hat = Panels[PanelID].GlobaltoLocal(hw.e_hat()); //photon polarization
	Vector3D n_hat = Panels[PanelID].Normal();

	// Get normal to plane of incidence (if n_hat and k_hat parallel, then ni_hat = y_hat)
	Vector3D ni_hat = abs(n_hat.dot(k_hat)) == 1 ? Vector3D(0, -1, 0) : n_hat.cross(k_hat);
	ni_hat.normalize();

	// (s) and (p) componenets of E-field (unitary)
	double es = e_hat.dot(ni_hat);
	double ep = sqrt(1 - es*es);

	// Get incidence angles
	double cosTi = abs(k_hat.Z());//cos(theta)
	double sinTi = sqrt(1 - cosTi*cosTi); //sin(theta)
	double cosPi = sinTi < 1E-10 ? 1 : k_hat.X() / sinTi; // cos(phi)
	double sinPi = sinTi < 1E-10 ? 0 : k_hat.Y() / sinTi; // sin(phi)

	cosTi = cosTi != 0.0 ? cosTi : 1E-15; // make cosTi a small number if zero
	sinTi = sinTi != 0.0 ? sinTi : 1E-15; // make sinTi a small number if zero
	cdouble kx1 = (n1*w/(double)SPEEDOFLIGHT) *sinTi;

	// Get transmition angles (complex)
	cdouble sinTt = n1*sinTi / n2; // sin(theta) // Snell's Law
	cdouble cosTt = sqrt(1.0 - sinTt*sinTt); // cos(theta)
	sinTt = sinTt != 0.0 ? sinTt : 1E-15; // make sinTt a small number if zero
	cosTt = cosTt != 0.0 ? cosTt : 1E-15; // make cosTt a small number if zero

	/*****************************  Fresnel Coefficients (s-polarization)**********************************/
	// get transfer matrix for multilayer film
	m11 = 1.0; m12 = 0.0; m21 = 0.0; m22 = 1.0;
	if (NumLayers > 0) TransferMatrix(kx1, w, kz_Dir, "s", m11, m12, m21, m22);

	// calculate the reflection and transmission
	p1 = - Z0*mu1 / (n1*cosTi);
	p2 = - Z0*mu1 / (n2*cosTt);
	ct = 1.0;
	cdouble rs = ((m11*p2 + m12) - (m21*p2 + m22)*p1) /
		         ((m11*p2 + m12) + (m21*p2 + m22)*p1);
	cdouble ts = 2.0*ct*p2 /
		         ((m11*p2 + m12) + (m21*p2 + m22)*p1);

	/*****************************  Fresnel Coefficients (p-polarization)**********************************/
	// get transfer matrix for multilayer film
	m11 = 1.0; m12 = 0.0; m21 = 0.0; m22 = 1.0;
	if (NumLayers > 0) TransferMatrix(kx1, w, kz_Dir, "p", m11, m12, m21, m22);
	
	// calculate the reflection and transmission
	p1 = Z0*mu1*cosTi/ n1;
	p2 = Z0*mu2*cosTt/ n2;
	ct = cosTi / cosTt;
	cdouble rp = ((m11*p2 + m12) - (m21*p2 + m22)*p1) /
		         ((m11*p2 + m12) + (m21*p2 + m22)*p1);
	cdouble tp = 2.0*ct*p2 /
		         ((m11*p2 + m12) + (m21*p2 + m22)*p1);

	//**************************** Effective reflectivity and transmissivity
	// reflected and transmitted power s-polarization
	double Rs = abs(rs*conj(rs));
	double Ts = real(n2*cosTt) / real(n1*cosTi)*abs(ts*conj(ts));
	// reflected and transmitted power p-polarization
	double Rp = abs(rp*conj(rp));
	double Tp = real(conj(n2)*cosTt) / real(conj(n1)*cosTi)*abs(tp*conj(tp));

	// Get effective reflectivity and transmisivity
	double R = Rs*abs(es*conj(es)) + Rp*abs(ep*conj(ep));
	double T = Ts*abs(es*conj(es)) + Tp*abs(ep*conj(ep));
	double A = 1 - R - T;

	/******************************Set new photon state **********************************************/
	// Decide whether the photon gets reflected or transmitted
	//srand(int(time(NULL))*omp_get_thread_num()); // generate a seed
	//bool newstate = false;
	int kz_sgn = k_hat.Z() > 0 ? 1 : -1; // +- direction of kz
	Vector3D new_k, new_e;

	double Rand = RandomNum;

	// Keep polarizations (respecto to k-vector)
	new_e.SetVector(ep, es);

	if (Rand <= R) {							// Photon reflected
		// Photon keeps same kx and ky and flips kz
		new_k.SetVector(k_hat.X(), k_hat.Y(), -k_hat.Z());
		kz_sgn *= -1; //change direction of kz
		FresnelState = 1;
	}
	else if (Rand > R && Rand <= R + T) {		// Photon Transmited
		// Photon keeps same kx and ky. kz direction is conserved, 
		// but changes magnitude based on snell's law
		new_k.SetVector(k_hat.X(), k_hat.Y(), real(cosTt)*kz_sgn);
		FresnelState = -1;
	}
	else if (Rand > R + T && Rand <= R + T + A) {	// Photon Absorbed
		new_k = k_hat;
		FresnelState = 0;
	}

	if (FresnelState != 0) {
		// *** Set new possition of photon (Local coordintes)
		// Move photon to the other side of the surface by a constant
		// (moves photon over 10xTolerance to avoid being detected by the same panel)
		double MoveConst = abs(new_k.Z())<Tol ? 10 : abs(10 * Tol / new_k.Z());
		x0_loc.x += MoveConst*new_k.X();
		x0_loc.y += MoveConst*new_k.Y();
		x0_loc.z += 10 * Tol*kz_sgn;

		// *** Store results
		// Store photon's new possition (pass to global coordinates) 
		hw.SetPosition(Panels[PanelID].LocaltoGlobal(x0_loc));

		// Store photon's momentum (pass to global coordinates)
		hw.SetMomentum(Panels[PanelID].LocaltoGlobal(new_k));

		// Rotate pol. to global coordinates and store new polarization
		hw.SetPolarization(new_e);
	}
	return FresnelState; // returns the final state of the photon
}


/**************************************** Subclass RTMonitor ****************************************/
RTMonitor::RTMonitor() : Surface() { 
	mon_initialize();
}

RTMonitor::RTMonitor(RTRegion *External, RTRegion *Internal, string &GeometryType, 
	const string &monLabel, const int Nhw, string plot_type, const int &idx) :
	Surface(External, Internal, GeometryType, monLabel, idx) {
	
	check_plot_type(plot_type);
	mon_initialize();
	TotalPhotons = Nhw;
}

RTMonitor::~RTMonitor() {
	Close_oFiles();
	mon_initialize();
	CleanUp();
}

void RTMonitor::mon_initialize(void) {
	foutID = NULL;
	x_data = NULL; x_size = 0;
	y_data = NULL; y_size = 0;
	print_to_file = true;

	//fextra = NULL;
}

void RTMonitor::Close_oFiles(void) {
	if (foutID != NULL) fclose(foutID);
	if (x_data != NULL) delete[] x_data;
	if (y_data != NULL) delete[] y_data;
	mon_initialize();

	//if (fextra != NULL) fclose(fextra);
}

bool RTMonitor::PhotonHits(Photon &hw) const {
	int iPanel;
	bool IsHit = false;
	Point3D x0, x0_new;
	double dz = DBL_MAX, dz_new;

		for (int i = 0; i < NumPanels; i++) {
			if (Panels[i].PhotonHits(hw, dz_new, x0_new)) {
				if (dz_new < dz) {
					IsHit = true;
					iPanel = i;
					dz = dz_new;
					x0 = x0_new;
				}
			}
		}
		if (IsHit)	hw.SetMonitorHit(x0, dz, iPanel, Index);
	return IsHit;
}

void RTMonitor::check_plot_type(string plot_type) {
	if (!plot_type.compare("RADIATION PROPERTIES")) out_type = RAD_PROP;
}

void RTMonitor::AddTotalPhoton(const int Nhw) {
	TotalPhotons = Nhw;
}

dvector* RTMonitor::Get_Output(void) const {
	return y_data;
}

int RTMonitor::GetPhotonCount(void) const {
	return TotalPhotons;
}

void RTMonitor::SetOutput(double *w_freq, int w_size,  string msh_file_name) {
/* Function set number of output results
	Definition:
		output 0: Radiation properties
		output 1: Fresnel coefficients
		output 2: Power distribution
*note: rest of output file will be defined in the future
*/
	string oFile;

	switch (out_type) {
	case RAD_PROP: // just radiation properties for now
		// create output data arrays
		y_size = 3;
		// if print to file is true, create file
		if (print_to_file) {
			oFile = Label + ".pow";
			foutID = fopen(oFile.c_str(), "w");

			fprintf(foutID, "# Radiation power file\n");
			fprintf(foutID, "# File columns: \n");
			fprintf(foutID, "# \t 1 Incidence zenith angle (deg) \n");
			fprintf(foutID, "# \t 2 Incidence azimuth angle (deg) \n");
			fprintf(foutID, "# \t 3 Frequency (rad/s)\n");
			fprintf(foutID, "# \t 4 Photon's relative power intensity \n");
			fprintf(foutID, "# \t 5 Photon's relative energy flux \n");
			fprintf(foutID, "# \t 6 Mean distance (um) traveled by photon \n");
			fprintf(foutID, "# \n");

			// extra output
			//fextra = fopen(Label.c_str(), "w");
		}
		break;
	default:
		ErrorMsg("Monitor: Unrecognized output");
		break;
	}
	set_x_data(w_freq, w_size, msh_file_name);
	set_y_data();
}

void RTMonitor::set_x_data(double *w_freq, int w_size, string msh_file_name) {
	switch (out_type) {
	case RAD_PROP:
		x_size = 1;
		x_data = new dvector[x_size];

		// save frequency data list
		for (int ix = 0; ix < w_size; ix++)
			x_data[0].push_back(w_freq[ix]);
		break;
	default:
		ErrorMsg("Monitor: Unrecognized output");
	}
}

void RTMonitor::set_y_data(void) {
	int simulation_size = 1;
	// Determine number of output values to be stored
	for (int ix = 0; ix < x_size; ix++)
		simulation_size *= (int)x_data[ix].size();

	// generate output arrays and set all values to 0
	y_data = new dvector[y_size];
	for (int iy = 0; iy < y_size; iy++)
		for (int ix_dat = 0; ix_dat < simulation_size; ix_dat++)
			y_data[iy].push_back(0.0);
}

void RTMonitor::clean_y_data(void) {

	for (int iy = 0; iy < y_size; iy++)
		for (int id = 0; id < y_data[iy].size(); id++)
			y_data[iy].at(id) = 0.0;
}

void RTMonitor::PhotonDetected(const int &w_idx,Photon *hw) {
	int iw;
	double Tol = 1E-4;
	int PanelID = hw->GetMonitorPanelID();
	Vector3D k_hat = Panels[PanelID].GlobaltoLocal(hw->k_hat()); // photon momentum in local
	Vector3D n_hat = Panels[PanelID].Normal(); // normal to monitor panel
	double cosTheta = abs(n_hat.dot(k_hat)); // Theta: angle of incidence of photon
	Point3D x0_loc = hw->x_MonitorHits(); // location where photon hits panel
	int kz_sgn = SIGN(k_hat.Z()); // +- direction of kz

	double MoveConst = abs(k_hat.Z()) < Tol ? 10 : abs(10 * Tol / k_hat.Z());
	x0_loc.x += MoveConst * k_hat.X();
	x0_loc.y += MoveConst * k_hat.Y();
	x0_loc.z += 10 * Tol*kz_sgn;

	// *** Store results
	// Store photon's new possition (pass to global coordinates) 
	hw->SetPosition(Panels[PanelID].LocaltoGlobal(x0_loc));

	// Photon is not counted if comes oposite to monitor's direction
	if (!hw->PhotonIs(Regions[0])) return; // Photon from Region[1] to Region[0] >> not detected!

	// Photon crosses from Region[0] to Region[1] >>> detected!!

	switch (out_type) {
	case RAD_PROP:
		/* If photon is emited from the source (not from a QD) the frequency index matches
		that from the source. Otherwise, a search function has to be activated to find the
		correding freq index from w_array*/
		
		if (!hw->is_QDemitted()) iw = w_idx; // get freq index form the source
		// else // "iw" from search in the w_array

		y_data[0].at(iw) += 1.0 / (double)TotalPhotons;
		y_data[1].at(iw) += cosTheta / (double)TotalPhotons;
		y_data[2].at(iw) += hw->get_traveled_path()/ (double)TotalPhotons;

		//fprintf(fextra, "%7.3e ", hw->get_traveled_path());
		break;
	default:
		ErrorMsg("Monitor: Unrecognized output");
	}
}


void RTMonitor::set_monitor(const double &xtheta, const double &xphi) {
	theta = xtheta;
	phi = xphi;
}
void RTMonitor::set_file_output_off(void) {
	print_to_file = false;
}

bool RTMonitor::Print_To_File(void) const {
	return print_to_file;
}

void RTMonitor::save_to_file(void) {
	//double P;
	if (foutID == NULL || !print_to_file) 
		ErrorMsg("Error Monitor " + Label + ": file not set");

	switch (out_type) {
	case RAD_PROP:
		for (int iw = 0; iw < x_data[0].size(); iw++) {
			double Detect = y_data[0].at(iw) == 0 ? 1E-5 : y_data[0].at(iw);
			fprintf(foutID, "%-6.2f\t", theta); // Zenith
			fprintf(foutID, "%-6.2f\t", phi); // Azimuth
			fprintf(foutID, "%-7.4e\t", x_data[0].at(iw)); // frequency
			fprintf(foutID, "%-9.6e\t", y_data[0].at(iw)); // Net power out
			fprintf(foutID, "%-9.6e\t", y_data[1].at(iw)/ Detect); // Energy flux out
			fprintf(foutID, "%10.6e\n",  y_data[2].at(iw)/ Detect); // Av. distance traveled
		}
	}
}

dvector* RTMonitor::get_x_data(int &x_size_out) const{
	x_size_out = x_size;
	return x_data;
}

dvector* RTMonitor::get_y_data(int &y_size_out) const {
	y_size_out = y_size;
	return y_data;
}

double* RTMonitor::y_data_at(const int &idx) const {
	double *y_out = new double[y_size];
	for (int i = 0; i < y_size; i++) {
		y_out[i] = y_data[i].at(idx);
	}
	return y_out;
}

/**************************************** Main Class Surface ****************************************/
Surface::Surface()
{
	Initialize();
}

Surface::Surface(RTRegion *External, RTRegion *Internal, string &GeometryArg, 
	const string &xLabel, const int &idx)
{
	Initialize();
	SetGeometry(GeometryArg, idx);
	Label = xLabel;
	Regions[0] = External;
	Regions[1] = Internal;
}

Surface::~Surface()
{
	CleanUp();
}

void Surface::Initialize()
{
	vtx = NULL;
	Panels = NULL;
	nVertex = 0;
	NumPanels = 0;
	MeshTag = -1; // Part of ReadGMSHFile. Read the entire *.msh file. To be modified in the future
}

void Surface::CleanUp()
{
	if (NULL != vtx) delete[] vtx;
	if (NULL != Panels) delete[] Panels;
	Initialize();
}

void Surface::SetGeometry(string &GeometryArg, const int &idx) {
	Index = idx;
	string MeshFile;
	Point3D cdisplace = center;
	double Lbox[3];
	string SurfaceID = "Surface (" + to_string(Index) + ")";

	objArgsearch GeometryInst[] = {
		{"MESHFILE"			, PA_STRING , false, 1, 1, (void *)&MeshFile	},
		{"BOX"				, PA_DOUBLE , false, 3, 3, (void *)&Lbox		},
		{"BOX_Z+OPEN"		, PA_DOUBLE , false, 3, 3, (void *)&Lbox		},
		{"BOX_Z-OPEN"		, PA_DOUBLE , false, 3, 3, (void *)&Lbox		},
		{"BOX_ZOPEN" 		, PA_DOUBLE , false, 3, 3, (void *)&Lbox		},
		{"FLATSURFACE_Z"	, PA_DOUBLE , false, 2, 2, (void *)&Lbox		},
		{"CENTER"			, PA_POINT3D, false, 1, 1, (void *)&cdisplace  },
		{"",0,0,0,0,0}
	};

	GetArguments(GeometryArg, GeometryInst, 2, SurfaceID);
	if      (GeometryInst[0].IsFound){
		MeshFile = MeshFile + ".msh";
		FILE *MeshFileID = fopen(MeshFile.c_str(), "r");
		ReadGMSHFile(MeshFileID, MeshFile);
	}
	else if (GeometryInst[1].IsFound) {
		MakeBox(Lbox[0], Lbox[1], Lbox[2]);
	}
	else if (GeometryInst[2].IsFound) {
		MakeBox(Lbox[0], Lbox[1], Lbox[2], "Z(+)_open");
	}
	else if (GeometryInst[3].IsFound) {
		MakeBox(Lbox[0], Lbox[1], Lbox[2], "Z(-)_open");
	}
	else if (GeometryInst[4].IsFound) {
		MakeBox(Lbox[0], Lbox[1], Lbox[2], "Z_open");
	}
	else if (GeometryInst[5].IsFound) {
		MakeFlatSurface(Lbox[0], Lbox[1]);
	}
	SetCenter(); // get center of surface

	// Check if CENTER statement has been declared and move surface accordingly
	int i_end = sizeof(GeometryInst) / sizeof(GeometryInst[0]) - 2;
	if (GeometryInst[i_end].IsFound)
		MoveSurface(cdisplace);

	//cout << Print();
	//cout << "end Surface" << endl;

}

string Surface::Print() {
	string output = "";
	output  += "****************** Description of Panels *****************\n";
	for (int i = 0; i < NumPanels; i++) {
		output  += Panels[i].Print() + "\n";
	}

	output += "************** Distance of panels from the center *********\n";
	output += "Center:" + center.Print() + "\n";
	for (int i = 0; i < NumPanels; i++) {
		Vector3D Vest(Panels[i].Centroid() - center);
		output += "	n_hat*[r(" + to_string(i) + ") - c] = " +
			to_string(Vest.dot(Panels[i].Normal())) + "\n";
	}
	return output;
}

void Surface::MakeBox(double Lx, double Ly, double Lz, string open) {
	nVertex = 8;
	vtx = new Point3D[nVertex]; int i;
	vtx[0].SetPoint3D(  0,  0,  0);
	vtx[1].SetPoint3D( Lx,  0,  0);
	vtx[2].SetPoint3D( Lx, Ly,  0);
	vtx[3].SetPoint3D(  0, Ly,  0);
	vtx[4].SetPoint3D(  0,  0, Lz);
	vtx[5].SetPoint3D( Lx,  0, Lz);
	vtx[6].SetPoint3D( Lx, Ly, Lz);
	vtx[7].SetPoint3D(  0, Ly, Lz);

	if (!open.compare("Z(+)_open")) {
		NumPanels = 10;
		Panels = new Panel[NumPanels];
		i = 0;
		Panels[i].SetCoord(vtx[0], vtx[2], vtx[1]); i++;
		Panels[i].SetCoord(vtx[0], vtx[3], vtx[2]); i++;
		Panels[i].SetCoord(vtx[1], vtx[4], vtx[0]); i++;
		Panels[i].SetCoord(vtx[5], vtx[4], vtx[1]); i++;
		Panels[i].SetCoord(vtx[2], vtx[3], vtx[7]); i++;
		Panels[i].SetCoord(vtx[2], vtx[7], vtx[6]); i++;
		Panels[i].SetCoord(vtx[0], vtx[4], vtx[3]); i++;
		Panels[i].SetCoord(vtx[4], vtx[7], vtx[3]); i++;
		Panels[i].SetCoord(vtx[1], vtx[2], vtx[5]); i++;
		Panels[i].SetCoord(vtx[5], vtx[2], vtx[6]); i++;
	}
	else if (!open.compare("Z(-)_open")) {
		NumPanels = 10;
		Panels = new Panel[NumPanels];
		i = 0;
		Panels[i].SetCoord(vtx[1], vtx[4], vtx[0]); i++;
		Panels[i].SetCoord(vtx[5], vtx[4], vtx[1]); i++;
		Panels[i].SetCoord(vtx[2], vtx[3], vtx[7]); i++;
		Panels[i].SetCoord(vtx[2], vtx[7], vtx[6]); i++;
		Panels[i].SetCoord(vtx[0], vtx[4], vtx[3]); i++;
		Panels[i].SetCoord(vtx[4], vtx[7], vtx[3]); i++;
		Panels[i].SetCoord(vtx[1], vtx[2], vtx[5]); i++;
		Panels[i].SetCoord(vtx[5], vtx[2], vtx[6]); i++;
		Panels[i].SetCoord(vtx[5], vtx[6], vtx[4]); i++;
		Panels[i].SetCoord(vtx[6], vtx[7], vtx[4]); i++;
	}
	else if (!open.compare("Z_open")) {
		NumPanels = 8;
		Panels = new Panel[NumPanels];
		i = 0;
		Panels[i].SetCoord(vtx[1], vtx[4], vtx[0]); i++;
		Panels[i].SetCoord(vtx[5], vtx[4], vtx[1]); i++;
		Panels[i].SetCoord(vtx[2], vtx[3], vtx[7]); i++;
		Panels[i].SetCoord(vtx[2], vtx[7], vtx[6]); i++;
		Panels[i].SetCoord(vtx[0], vtx[4], vtx[3]); i++;
		Panels[i].SetCoord(vtx[4], vtx[7], vtx[3]); i++;
		Panels[i].SetCoord(vtx[1], vtx[2], vtx[5]); i++;
		Panels[i].SetCoord(vtx[5], vtx[2], vtx[6]); i++;
	}
	else {
		NumPanels = 12;
		Panels = new Panel[NumPanels];
		i = 0;
		Panels[i].SetCoord(vtx[0], vtx[2], vtx[1]); i++;
		Panels[i].SetCoord(vtx[0], vtx[3], vtx[2]); i++;
		Panels[i].SetCoord(vtx[1], vtx[4], vtx[0]); i++;
		Panels[i].SetCoord(vtx[5], vtx[4], vtx[1]); i++;
		Panels[i].SetCoord(vtx[2], vtx[3], vtx[7]); i++;
		Panels[i].SetCoord(vtx[2], vtx[7], vtx[6]); i++;
		Panels[i].SetCoord(vtx[0], vtx[4], vtx[3]); i++;
		Panels[i].SetCoord(vtx[4], vtx[7], vtx[3]); i++;
		Panels[i].SetCoord(vtx[1], vtx[2], vtx[5]); i++;
		Panels[i].SetCoord(vtx[5], vtx[2], vtx[6]); i++;
		Panels[i].SetCoord(vtx[5], vtx[6], vtx[4]); i++;
		Panels[i].SetCoord(vtx[6], vtx[7], vtx[4]); i++;

	}
}

// function only works on z-oriented surfaces. Needs to be improved for to get any surface orientation
void Surface::MakeFlatSurface(double Lx, double Ly, string n_hat) {
	nVertex = 4;
	vtx = new Point3D[nVertex]; int i;

	NumPanels = 2;
	Panels = new Panel[NumPanels];
	
	if (!n_hat.compare("z")){
		vtx[0].SetPoint3D(0, 0, 0);
		vtx[1].SetPoint3D(Lx, 0, 0);
		vtx[2].SetPoint3D(Lx, Ly, 0);
		vtx[3].SetPoint3D(0, Ly, 0);
	}

	i = 0;
	Panels[i].SetCoord(vtx[0], vtx[2], vtx[1]); i++;
	Panels[i].SetCoord(vtx[0], vtx[3], vtx[2]); i++;
}

void Surface::SetCenter(void) {
	if (vtx != NULL){
		Point3D newcenter(0, 0, 0);
		for (int i = 0; i < nVertex; i++) {
			newcenter += vtx[i];
		}
		center = newcenter / nVertex;
	}
	else {
		ErrorMsg("Error: Undefined vertex for surface: " + Label);
	}
}

void Surface::MoveSurface(Point3D r) {
	for (int i = 0; i < NumPanels; i++)
		Panels[i].MovePanel(r - center);
	for (int i = 0; i < nVertex; i++)
		vtx[i] = r - vtx[i];
	center = r;
}

int Surface::GetNumPanels(void) const {
	return NumPanels;
}

Panel* Surface::GetPanel(int i) const {
	return &Panels[i];
}
