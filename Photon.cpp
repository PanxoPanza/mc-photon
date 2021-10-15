#define _USE_MATH_DEFINES
#include "mcphoton_lib.h"
#define EXACT_DIST 0
#define HG_APPROX  1

Photon::Photon()
{
	Initialize();
	Allocate();
}

Photon::Photon(double w)
{
	Initialize();
	Allocate();
	frequency = w;
}

Photon::Photon(double w, Point3D x)
{
	Initialize();
	Allocate();
	frequency = w;
	SetPosition(x);
}

Photon::Photon(double w, Point3D x, Vector3D M)
{
	Initialize();
	Allocate();
	frequency = w;
	SetPosition(x);
	SetMomentum(M);
}

Photon::Photon(double w, Point3D x, Vector3D M, Vector3D P)
{
	Initialize();
	Allocate();
	frequency = w;
	SetPosition(x);
	SetMomentum(M);
	SetPolarization(P);
}

Photon::~Photon(void) { CleanUp(); }

void Photon::Initialize(void)
{
	Position = NULL;
	Momentum = NULL;
	Polarization = NULL;
	PhotonRegion = NULL;
	frequency = 0;
	Alive = 1;
	QYLoss = 0;
	NumReflects = 0;
	oldRegionLabel = "";
	QD_emitted = false;
	IsExterior = true; // photon is assumed to be in the exterior by default
	TravelPath = 0;

	// Initialize path length and other temporary parameters
	i_particle = -1;

	dx_Boundary = DBL_MAX / 2;

	dx_Surface = DBL_MAX;
	SurfaceID = -1;
	SurfacePanelID = -1;
	SurfacePanelID = -1;

	dx_Monitor = DBL_MAX;
	MonitorPanelID = -1;

	ExtLength = DBL_MAX;
	ParticleLength = DBL_MAX;
}

void Photon::InitializeLength(void) {

	i_particle = -1;

	dx_Boundary = DBL_MAX/2;

	dx_Surface = DBL_MAX;
	SurfaceID = -1;
	SurfacePanelID = -1;
	SurfacePanelID = -1;

	dx_Monitor = DBL_MAX;
	MonitorPanelID = -1;

	// set properties related with the Region
	if (PhotonRegion == NULL) // check if region has been defined
		ErrorMsg("Undefined Region in photon object");

	SetParticlePathlength(); // set path length to a particle
	SetExtLength();			 // set path length to absorption
	set_travel_path();       // set distance traveled by the photon on previous event
}

void Photon::set_travel_path(void) {
	if (Position == NULL) ErrorMsg("Undefined photon position");

	/* if previous region was not exterior (IsExterior = false) track the traveled distance */
	if (!IsExterior) TravelPath += x_old.distance(*Position);

	// Photon is in exterior (true), else (false)
	IsExterior = (bool)!PhotonRegion->Label.compare("EXTERIOR");
	//printf("TravelPath: %8.3f, Region: %s\n", TravelPath, PhotonRegion->Label.c_str());
	//cout << Position->Print() << endl;
	x_old = *Position;
}

void Photon::CleanUp(void)
{
	if (Position != NULL) delete Position;
	if (Momentum != NULL) delete Momentum;
	if (Polarization != NULL) delete Polarization;
	if (!fwTheta.isempty()) fwTheta.clear();
	// Don't do "delete PhotonRegion";  would delete Region's object!!
	Initialize();
}

void Photon::Allocate(void)
{
	Position = new Point3D;
	Momentum = new Vector3D;
	Polarization = new Vector3D;
}

void Photon::SetRegion(RTRegion *FromRegion) {
	if (PhotonRegion == NULL) {// First time region designation
		PhotonRegion = FromRegion;
		IsExterior = (bool)!PhotonRegion->Label.compare("EXTERIOR");
		if (Position == NULL) ErrorMsg("Undefined photon's position");
		x_old = *Position;
	}
	else // if a new region is set
		if (PhotonRegion->Label.compare(FromRegion->Label) != 0)
			PhotonRegion = FromRegion;
}

void //Sets Photon momentum vector
Photon::SetMomentum(Vector3D M){
    *Momentum = M.normalize();
}

void //Sets Photon Position vector
Photon::SetPosition(Point3D X){
    *Position = X;
}

void //SetPolarization Vector
Photon::SetPolarization(Vector3D P){
	P = k_LocaltoGlobal(P);
    *Polarization = P.normalize();
}

void //Set if photon has been killed by QY
Photon::SetQYLoss() {
	QYLoss = 1;
}

void Photon::SetSurfaceHit(const Point3D &x,const double &D,
	                       const int &iPanel, const int &iSurface) {
	if (D < dx_Surface){
		x_hitSurface = x;
		dx_Surface = D;
		SurfaceID = iSurface;
		SurfacePanelID = iPanel;
	}
}

void Photon::SetMonitorHit(const Point3D &x, const double &D,
							const  int &iPanel, const int &iMonitor) {
	if (D < dx_Monitor) {
		x_hitMonitor = x;
		dx_Monitor = D;
		MonitorID = iMonitor;
		MonitorPanelID = iPanel;
	}
}

void Photon::SetParticlePathlength(void) {

	if (PhotonRegion != NULL){ // If the region has been defined
		if (PhotonRegion->IsComposite()) {// if the material has nanoparticles
			ParticleLength = PhotonRegion->NP_sample_path(i_particle);
			fwTheta = PhotonRegion->NP_wfDist_at(i_particle);
			switch (PhotonRegion->NP_PhaseDistType_at(i_particle)) {
			case EXACT_DIST: // Exact distribution from file
				fwTheta = PhotonRegion->NP_wfDist_at(i_particle);
				break;

			case HG_APPROX: // Henyney-Greenstein approximation
				gasym = PhotonRegion->NP_wAsymmetry_at(i_particle);
				break;

			default:
				ErrorMsg("Unknown phase distribution function in Region: " + PhotonRegion->Label);
			}
		}
		else
			ParticleLength = DBL_MAX;
	}
	else {
		cout << "Photon, Error! Region has not been set" << endl;
		exit(EXIT_FAILURE);
	}
}

void Photon::SetExtLength(void) {

	if (PhotonRegion != NULL) {
		double MeanPathExt = PhotonRegion->get_wExtPath();
		if (MeanPathExt < DBL_MAX) ExtLength = -MeanPathExt * log(RandomNum);
		else ExtLength = DBL_MAX;
	}
	else {
		cout << "Photon, Error! Region has not been set" << endl;
		exit(EXIT_FAILURE);
	}
}

// Photon is scattered(true) or absorbed(false) by a particle
bool Photon::HitsParticle(void) {
	if (RandomNum <= PhotonRegion->NP_wScaAlbedo_at(i_particle)) {
		NPscatter(); return true;
	}

	NPabsorb();
	return false;
}

double Photon::SpinTheta(void)
{
	double costheta, mu;
	switch(PhotonRegion->NP_PhaseDistType_at(i_particle)) {

		case EXACT_DIST : // Exact distribution from file
			costheta = fwTheta(RandomNum);
			break;

		case HG_APPROX  : // Henyney-Greenstein approximation
			if (gasym == 0) costheta = 2*RandomNum - 1;  // isotropic
			//else if (gasym == 1) costheta = 1;
			else {
				mu = (1 - gasym * gasym) / (1 - gasym + 2.0*gasym*RandomNum);
				costheta = (1 + gasym * gasym - mu * mu) /  (2.0 *gasym);
			  if    (costheta < -1) costheta = -1;
			  else if(costheta > 1) costheta = +1;
			}
			break;

		default :
			ErrorMsg("Unknown phase distribution function in Region: " + PhotonRegion->Label);
	}
	return costheta;
}

void Photon::NPscatter(void) {
	Vector3D nz_hat, ns_hat, new_e;
	Point3D new_pos;
	double u, v, w;
	double u1, v1, w1, temp, mu;
	double costheta, sintheta;
	double cosphi, sinphi;
	double Vol = PhotonRegion->NP_volume_at(i_particle); // Volume of particle (um)
	double Dp = 2 * pow(3 * Vol / (4 * M_PI), 1 / 3.); // effective diameter of the particle (um)

	u = Momentum->X();
	v = Momentum->Y();
	w = Momentum->Z();

	// get possition of photon at scatter point
	new_pos = *Position;
	new_pos.x += ParticleLength *u;
	new_pos.y += ParticleLength *v;
	new_pos.z += ParticleLength *w;

	// Get costheta from phase distribution function
  costheta = SpinTheta();
	sintheta = sqrt(1.0 - costheta * costheta);

	/* Sample psi. */
	double phi = 2.0*M_PI*RandomNum;
	cosphi = cos(phi);
	if (phi < M_PI)
		sinphi = + sqrt(1.0 - cosphi * cosphi);
	else
		sinphi = - sqrt(1.0 - cosphi * cosphi);

	/* New trajectory. */
	if (fabs(w) > 1 - 1E-12) {
		u1 = sintheta * cosphi;
		v1 = sintheta * sinphi;
		w1 = costheta*SIGN(w);
	}
	else {
		temp = sqrt(1 - w * w);
		u1 = costheta * u + sintheta * (cosphi*u*w - sinphi*v)/temp;
		v1 = costheta * v + sintheta * (cosphi*v*w + sinphi*u)/temp;
		w1 = costheta * w - sintheta*cosphi*temp;
	}
	u = u1;
	v = v1;
	w = w1;

	// normal vector pointing to z direction
	///nz_hat.SetVector(0, 0, 1);

	// Get normal of the scattering plane
	//ns_hat = abs(nz_hat.dot(*Momentum)) == 1 ?
	//					Vector3D(0, -1, 0) : nz_hat.cross(*Momentum);
	//ns_hat.normalize(); // normalize ns_hat

	// get E-components respect to the scattering plane
	//double es = Polarization->dot(ns_hat);
	double es = RandomNum;
	double ep = sqrt(1 - es*es);
	new_e.SetVector(ep, es, 0); 		// save intial polarization

	// save new values (update postion outside of the particle's surface)
	SetPosition(new_pos); 					// new position
	SetMomentum(Vector3D(u, v, w)); // new momentum
	SetPolarization(new_e); 				// new polarization
}

void Photon::NPabsorb(void) { // absorb the photon and kill it
	PhotonKill();
}

Vector3D //Get Momentum Vector
Photon::k_hat() const {
    return *Momentum;
}

Point3D //Get Position
Photon::x0() const {
	return *Position;
}

RTRegion* Photon::GetRegion() const {
	return PhotonRegion;
}

bool Photon::PhotonIs(const RTRegion *Reg) const{
	string RegionLabel = PhotonRegion->Label;
	if (RegionLabel.compare(Reg->Label) == 0) return true;
	return false;
}

// Rotates a unitary vector based on momentum's local coordinates
// into a global coordinate system
Vector3D Photon::k_LocaltoGlobal(const Vector3D &V)
{
	double cosT = Momentum->Z();
	double sinT = sqrt(1 - cosT*cosT);
	double cosP = sinT < 1E-10 ? 1 : Momentum->X() / sinT;
	double sinP = sinT < 1E-10 ? 0 : Momentum->Y() / sinT;

	//Rotation of
	Vector3D U(
		+cosP*cosT*V.X() - sinP*V.Y() + cosP*sinT*V.Z(),
		+sinP*cosT*V.X() + cosP*V.Y() + sinP*sinT*V.Z(),
		     -sinT*V.X() +  0 * V.Y() +      cosT*V.Z());
	return U.normalize(); // returns a unitary vector
}

std::string Photon::Print(void) const {
	ostringstream Freq, Length[6];
	Freq << frequency;
	Length[0] << dx_Boundary;
	Length[1] << dx_Surface;
	Length[2] << dx_Monitor;
	Length[3] << ExtLength;
	Length[4] << ParticleLength;
	Length[5] << TravelPath;

	string outstring = "";
	outstring += "	Frequency =	" + Freq.str() + " rad/s\n";
	outstring += "	Position =	" + Position->Print() + "\n";
	outstring += "	Momentum =	" + Momentum->Print() + "\n";
	outstring += "	Polarization =	" + Polarization->Print() + "\n";
	outstring += "	Region =	" + PhotonRegion->Label + "\n";
	outstring += "	PathLengths (um):\n";
	outstring += "		Boundary	=	" + Length[0].str() + "\n";
	outstring += "		Surface		=	" + Length[1].str() + "\n";
	outstring += "		Monitor		=	" + Length[2].str() + "\n";
	outstring += "		Extinction	=	" + Length[3].str() + "\n";
	outstring += "		Particle	=	" + Length[4].str() + "\n";
	outstring += "		Travel Dist.	=	" + Length[5].str() + "\n";
	outstring += "		R. Count	=	" + to_string(NumReflects) + "\n";
	outstring += "	IsAlive =	" + to_string(Alive) + "\n";
	return outstring;
}

bool Photon::IsAlive() const{
    return Alive;
}

void //Kill Photon
Photon::PhotonKill(){
    Alive = 0;
}

bool //Get if photon has been killed by QY
Photon::GetQYLoss() const {
	return QYLoss;
}

Vector3D //Get Polarization Vector
Photon::e_hat() const {
	return *Polarization;
}

double //Get Wavelength
Photon::GetFrequency() const {
	return frequency;
}

double Photon::GetParticleLength(void) const { return ParticleLength; }
double Photon::GetExtLength() const {return ExtLength; }
string Photon::GetParticleLabel() const { return PhotonRegion->NP_label_at(i_particle); }

double Photon::GetBoundaryDistace(void) const { return dx_Boundary; }

double Photon::GetSurfaceDistance(void) const { return dx_Surface; }
int Photon::GetSurfaceID(void) const { return SurfaceID; }
int Photon::GetSurfacePanelID(void) const { return SurfacePanelID; }
Point3D Photon::x_SurfaceHits(void) const { return x_hitSurface;  }

double Photon::GetMonitorDistance(void) const { return dx_Monitor; }
int Photon::GetMonitorID(void) const { return MonitorID; }
int Photon::GetMonitorPanelID(void) const { return MonitorPanelID; }
Point3D Photon::x_MonitorHits(void) const { return x_hitMonitor; }
bool Photon::is_QDemitted(void) const { return QD_emitted; }

double Photon::get_traveled_path(void) const { return TravelPath; }
