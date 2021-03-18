#ifndef __MC_RayTracing__MCRT_Library__
#define __MC_RayTracing__MCRT_Library__

#include<iostream>
#include<string>
#include<math.h>
#include<fstream>
#include<sstream>
#include<omp.h>
#include<chrono>
#include<ctime>
#include <time.h>
#include "materials.h"
#include "localmathlib.h"

#include <cfloat>
#include <cstring>
#include <cstdio>
#include <cstdlib>

#define SPEEDOFLIGHT 299792458 // m/s
#define EV 1.602176565e-19 // eV in J
#define EPS0 8.85418782e-12 // F/m
#define MU0 M_PI*4e-7
#define Z0 MU0*SPEEDOFLIGHT // Impedance of freespace
#define H_PLANK 6.62606957e-34
#define HBAR H_PLANK / (2.*M_PI)
#define RandomNum (double) RandomGen(1, 0, NULL)
#define SIGN(x) ((x)>=0 ? 1:-1)
#define InitRandomGen (double) RandomGen(0, 1.0, NULL)
#define MAX_REFLECTION 1e4

#define PA_DOUBLE  0
#define PA_INT     1
#define PA_STRING  2
#define PA_CDOUBLE 3
#define PA_POINT3D 4
#define PA_DBLSTR  5
#define PA_BOOL	   6

#define RAD_PROPERTIES 0

#define	OS_UNI 1
#define	OS_WIN 2


#ifdef __linux__
//linux code goes here
#define OS_SYSTEM 1
#elif _WIN32
// windows code goes here
#define OS_SYSTEM 2
#else

#endif


const cdouble II(0.0, 1.0);

using namespace std;

typedef vector<double> dvector;

typedef struct doublestr {
	double Val;
	string strfunc;
} doublestr;

typedef struct objArgsearch {
	string Name;
	int Type;
	bool IsFound;
	int nArgs_min;
	int nArgs;
	void *Storage;
} objArgsearch;


typedef struct InstrSet {
	string Name;
	int nArgs;
	int nArgsMax;
	string Value;
} InstrSet;

class Point3D
{
public:
	double x, y, z;
	Point3D();
	Point3D(double x1, double y1, double z1 = 0);

	Point3D operator+(const Point3D &coord);    //addition
	Point3D &operator+=(const Point3D &coord);  ////assigning new result to the vector
	Point3D operator-(const Point3D &coord);    //substraction
	Point3D &operator-=(const Point3D &coord);  //assigning new result to the vector
	Point3D &operator *=(double t);
	Point3D &operator /(double t);
	void SetPoint3D(double x1, double y1, double z1 = 0);
	double distance(const Point3D &coord);     // distance to "coord"
	string Print(void);    //display value of vectors
};

class Vector3D
{
private:
	double x, y, z;
public:

	//constructor
	Vector3D();
	Vector3D(double x1, double y1, double z1 = 0);
	Vector3D(Point3D r);
	Vector3D(const Vector3D &vec);    //copy constructor

	//Set if not constructed
	void SetVector(double x1, double y1, double z1 = 0);

	// vector operator
	Vector3D operator+(const Vector3D &vec);    //addition
	Vector3D &operator+=(const Vector3D &vec);  ////assigning new result to the vector
	Vector3D operator-(const Vector3D &vec);    //substraction
	Vector3D &operator-=(const Vector3D &vec);  //assigning new result to the vector
	Vector3D operator*(double value);    //multiplication
	Vector3D &operator*=(double value);  //assigning new result to the vector.
	Vector3D operator/(double value);    //division
	Vector3D &operator/=(double value);  //assigning new result to the vector
	Vector3D &operator=(const Vector3D &vec);
	double dot(const Vector3D &vec); //scalar dot_product
	Vector3D cross(const Vector3D &vec);    //cross_product
	double magnitude();     //magnitude of the vector
	Vector3D normalize();   //normalized vector
	double square(); //gives square of the vector

	//Display
	double distance(const Vector3D &vec);    //gives distance between two vectors
	double X() const; //return x
	double Y() const; //return y
	double Z() const; //return z
	string Print(void);    //display value of vectors
};

// sub class of RTRegion, with nanoparticles properties
class Mie_Inclusion {
private:
	double fV; // particle's volume fraction
	void Initialize();
	int NumFreq; // Number of frequencies for mie-data

	// code gets scattering properties
	bool Mie_calculate;

	// mie-scatt data
	bool MieData;		 // True if mie-scat data is available
	dataInterpol sCross; // (interpolation function) Scattering Cross section (um^2)
	dataInterpol aCross; // (interpolation function) Absorption Cross section (um^2)
	dataInterpol gAsymm; // (interpolation function) Asymmetry factor

	// effective properties
	bool eff_opt;		  // true if effective optical properties is available
	dataInterpol eps[2];	  // (interpolation function) effective eps
	dataInterpol mu[2];     // (interpolation function) effective mu

	// phase distribution data
	bool is_pDist_exat;   // True if exact phase distribution is available
	int Ntheta;			  // Number of angles for phase distribution
	dvector cosT;		  // cos(theta) array for phase distribution
	dataInterpol *fTheta; // (interpolation function) phase function distribution

	// ***********  interpolated variables
	bool frequency_set; // true if the properties have been set at a frequency "w"
	double wMeanPathParticle;
	double wScatAlbedo;
	double wgAsym;
	cdouble wEps;
	cdouble wMu;
	dataInterpol wNPftheta;

	// Internal function
	void Open_MieFile(string MieFileName[]); // extract data from *.mie file
	int Read_Parameters(ifstream &File); // read mie-parameters from *.mie file
	int Read_Data(ifstream &File); // read mie-data from *.mie file
	int setNPfTheta(const double &w); // set interpolation function for theta at frequency "w"

public:
	bool is_scatter; // true if scattering properties are considered (false for effective media)
	Mie_Inclusion(void);
	Mie_Inclusion(const string &str_inst, const string &RegionLabel);
	~Mie_Inclusion(void);
	void CleanUp(); // There is a memory leak when this function is kept private!!

	void set_object(const string &str_inst, const string &RegionLabel);

	// Particle parameters (minimum input)
	string RegionLabel;
	string NP_Label;	  // Particle's label
	string Host_Label;	  // Host material label
	double NPcon;	      // Particle's concentration (1/um^3)
	double NP_Vol;        // Particle's volume (um^3)
	double NP_Leff;		  // Particle's effective length (um)
	string MieFile;		  // Name of scattering data file
	int Phase_dist_type;  // type of prediction for phase distribution (from setup file)
	bool MieFileOpen;	  // true if Mie Scattering Data is open

	// control function for region
	int set_at_frequency(const double &w);
	void reset_wProperties(void);
	bool is_opt_properties(void) const;
	bool is_mie_data(void) const;
	bool is_phasedist(void) const;

	bool wProperties_set(void) const;
	double get_wMeanPath(void) const;
	double get_wScatAlbedo(void) const;
	double get_wAsym(void)const;
	cdouble get_wEps(void) const;
	cdouble get_wMu(void) const;
	dataInterpol get_wfDist(void) const;
};

class RTRegion {
private :
	//int i_np;						// Index of inclusions with minimum path length
	int Num_Inclusions;				// number of inclusion
	Mie_Inclusion *mie_particle;	// particle inclusions
	Material *Properties;			// Region's optical properties
	bool is_composite;
	bool frequency_set;    // boolean for properties set at a frequency
	double frequency;      // the frequency at which the properties have been extracted

	// ***********  interpolated variables
	cdouble wEps;
	cdouble wMu;
	double wExtPath;

	void Initialize();
	void CleanUp();
	void Allocate();

public:
	string Label;
	string OpMaterialName;

	// Contructor and destructor
	RTRegion();
	RTRegion(const string &RegionName, const string &MaterialName);
	~RTRegion();

	// Set region's material properties and label
	void SetRegion(const string &RegionName, const string &MaterialName);

	// Set nanoparticle's scattering and absorption data
	void Set_Inclusions(string strInst);

	// Get private values
	bool IsComposite(void) const; // true if material has nanoparticles
	cdouble get_wEps(void); // returns the dielectric constant at the set frequency
	cdouble get_wMu(void);  // returns the relative permitivity at the set frequency
	double get_wExtPath(void);

	double NP_sample_path(int &i_np);

	// get specific values of inclusions in region (get by index i);
	int get_NumInclusions(void) const;
	double NP_concentration_at(const int &i = 0) const; //Get nanoparticle's concentration
	double NP_volfrac_at(const int &i = 0) const; //Get nanoparticle's volume fraction
	double NP_volume_at(const int &i = 0) const;
	double NP_EffLength_at(const int &i = 0) const;
	string NP_label_at(const int &i = 0) const;

	// function for properties set at frquency "w"
	int setProperties(const double &w); // set RTRegion properties at "w"
	int NP_PhaseDistType_at(const int &i = 0) const; // type of phase distribution
	double NP_wMeanPath_at(const int &i = 0) const; //Particle's mean path
	double NP_wScaAlbedo_at(const int &i = 0) const; // Particle's scattering albedo
	double NP_wAsymmetry_at(const int &i = 0) const; // Asymmetry factor at freq w (rad/s)
	cdouble NP_wEps_at(const int &i = 0) const; //particles dielectric constant
	cdouble NP_wMu_at(const int &i = 0) const; //particles permeability
	dataInterpol NP_wfDist_at(const int &i = 0) const;
	void reset_wProperties(void);

	string Print(void) const; // Print relevant object data
};

class Photon {

private:
	// ************** Photon properties *******************************************************************
	int i_particle;   // index of particle hit by the photon
	double frequency; // Frequency in (rad/s)
	Point3D *Position;
	Vector3D *Momentum;
	Vector3D *Polarization;
	RTRegion *PhotonRegion;
	bool Alive; // Whether photon is keep in the system
	bool QYLoss;
	bool QD_emitted; // if photon has been emitted from a quantum dot (false by default)

	Vector3D k_LocaltoGlobal(const Vector3D &V);

	//****************** Particle and Material Interactions ***********************************************
	double ParticleLength; // Particle's path length (mm)
	double ExtLength; // Extinction path length for absorbing materials (mm)
	double gasym;
	dataInterpol fwTheta;
	double SpinTheta(void); // Spin zenith angle (costheta)
	void NPscatter(void);  // scatter photon
	void NPabsorb(void);   // absorb photon (kill)

	// ************* Monitor traveled path of photon inside the material ************************************
	double TravelPath; // path traveled by photon inside materials
	bool IsExterior; // true if photon is in the exterior
	Point3D x_old; // photon's previous possition from last event
	void set_travel_path(void); // calculate the distance traveled by the photon ****************************

	//*********************** Interface check parameters ****************************************************
	double dx_Boundary; // Distance to the domain's boundary (mm)

	double dx_Surface; // Distance to closest physical interface (mm)
	int SurfaceID, SurfacePanelID; // index of closest physical interface and panel
	Point3D x_hitSurface; // point in the surface's panel where photon hits

	double dx_Monitor; // Distance to closest monitor interface (mm)
	int MonitorID, MonitorPanelID;  // index of closest monitor interface and panel
	Point3D x_hitMonitor; // point in the monitor's panel where photon hits

	// ************************ Allocation and clean up variables *******************************************
	void Initialize(void);
	void CleanUp(void);
	void Allocate(void);

public:
	int NumReflects;
	string oldRegionLabel;

	Photon();
	// Constructor without frequecy units declaration (rad/s by default)
	Photon(double w);
	Photon(double w, Point3D x);
	Photon(double w, Point3D x, Vector3D M);
	Photon(double w, Point3D x, Vector3D M, Vector3D P);
	// Destructor
	~Photon();

	// ******************************  Special functions *************************************
	//void RandomMomentum();
	bool IsAlive() const;
	void PhotonKill();
	bool PhotonIs(const RTRegion *Reg) const;
	string Print(void) const;
	void InitializeLength(void);
	bool HitsParticle(void); // Photon is scattered(true) or absorbed(false) by a particle
	//double RandomGen(char Type, long Seed, long *Status);
	//double PhotonRand(void);

	// ***************************** Functions to set private values ***********************
	void SetRegion(RTRegion *FromRegion);
	void SetPosition(Point3D x);
	void SetMomentum(Vector3D M);
	void SetPolarization(Vector3D P);
	void SetQYLoss();
	void SetParticlePathlength(void);
	void SetExtLength(void); // Internal function. Sets the extincion path length

	// ***************************** Functions to show private values ************************
	Point3D x0() const;
	Vector3D k_hat() const;
	Vector3D e_hat() const;
	RTRegion* GetRegion() const;
	double GetFrequency() const;
	bool GetQYLoss() const;
	double GetParticleLength() const;
	double GetExtLength() const;
	string GetParticleLabel(void) const;
	bool is_QDemitted(void) const;
	double get_traveled_path(void) const;

	// *************************  Controllers of photon interface check ********************
	//	Boundary checks
	double GetBoundaryDistace() const;

	//	RTSurface checks
	void SetSurfaceHit(const Point3D &x, const double &D,
		               const  int &iPanel, const int &iSurface);
	double GetSurfaceDistance(void) const;
	int GetSurfaceID() const;
	int GetSurfacePanelID() const;
	Point3D x_SurfaceHits() const;

	//	RTMonitor checks
	void SetMonitorHit(const Point3D &x, const double &D,
		const  int &iPanel, const int &iSurface);
	double GetMonitorDistance() const;
	int GetMonitorID() const;
	int GetMonitorPanelID() const;
	Point3D x_MonitorHits() const;
};

class Panel
{
protected:
	double A_panel;
	Point3D *vtx, *vtx2D, *centroid;
	Vector3D *normal;
	double theta, phi;
	double Tol; // Error tolerance
	void SetPlaneNormal(void); // Calculates the plane's normal unitary vector
	void SetNormalAngles(void); // Calculates the plane's normal angles
	void SetCentroid(void);

	// Allocate and clean up variables
	void Initialize(void);
	void CleanUp(void);
	void Allocate(void);

public:
	int nVtx;
	int Index;
	Panel();
	Panel(Point3D r1, Point3D r2, Point3D r3);
	~Panel();

	// Special functions
	Point3D  GlobaltoLocal(Point3D rx) const;
	Vector3D GlobaltoLocal(Vector3D Vx) const;
	Point3D  LocaltoGlobal(Point3D rx) const;
	Vector3D LocaltoGlobal(Vector3D Vx) const;
	double ElementArea(Point3D r1, Point3D r2, Point3D r3) const;
	bool IsInside(const Point3D &r) const;
	bool PhotonHits(const Photon &hw, double &dXif, Point3D &Xf) const;
	double DistanceToPhoton(const Photon &hw) const;
	string Print(void) const;
	void MovePanel(const Point3D &r);

	//Set values
	void SetCoord(Point3D r1, Point3D r2, Point3D r3);

	//Show private values
	Point3D Vertex(int i) const;
	double NormalTheta(void) const;
	double NormalPhi(void) const;
	Point3D Centroid(void) const;
	Vector3D Normal(void) const;
};

class Surface
{
protected:
	int NumPanels;
	int nVertex;
	Point3D *vtx;
	Panel *Panels;
	int MeshTag;
	Point3D center;
	int NumRegions;
	RTRegion *Regions[2];
	// add vertex of the polygon;

	void Initialize();
	void CleanUp();
	void MakeBox(double Lx, double Ly, double Lz, string open="");
	void MakeFlatSurface(double Lx, double Ly, string n_hat = "z");
	void SetCenter(void);
	void MoveSurface(Point3D r);

public:
	string Label;
	int Index;
	Surface();
	Surface(RTRegion *External, RTRegion *Internal, string &GeometryArg,
		const string &xLabel="", const int &idx = 0);
	~Surface();
	void SetGeometry(string &GeometryType, const int &idx);
	void ReadGMSHFile(FILE *MeshFile, const string &FileName);
	int GetNumPanels(void) const;
	Panel* GetPanel(int i) const;
	string Print();

};

class RTSurface : public Surface {
private :
	int NumLayers;
	Material *matFilm;
	double *tFilm;
	cdouble *epsFilm;
	cdouble *muFilm;
	bool frequency_set;

public:

	RTSurface();
	RTSurface(RTRegion *External, RTRegion *Internal, string &GeometryArg,
		const string &xLabel, const int &idx = 0);
	~RTSurface();
	bool PhotonHits(Photon &hw) const;
	int RefractPhoton(Photon &hw) const;
	void SetCoating(string strThinFilms); // set multilayer coating
	void TransferMatrix(const cdouble &kxi, double w, int kz_Dir, string pol,
		cdouble &m11, cdouble &m12, cdouble &m21, cdouble &m22) const;
	int FresnelRefraction(const int &idx_from, const int &idx_to, const int &kz_Dir, Photon &hw) const;
	int set_wLayers(const double &w); // set optical properties of the layer at frequency "w"
	void reset_wProperties(void);
};

class RTMonitor : public Surface {
private :
	dvector *x_data; int x_size;
	dvector *y_data; int y_size;
	FILE *foutID; // , *fextra;
	int TotalPhotons;
	int out_type;
	bool print_to_file;
	void mon_initialize(void);
	void set_y_data(void);
	void set_x_data(double *w_freq = NULL, int w_size = 0, string msh_file_name = "");
	void check_plot_type(string plot_type);

public:
	double theta, phi;
	RTMonitor();
	RTMonitor(RTRegion *External, RTRegion *Internal, string &GeometryType,
		const string &monLabel ,const int Nhw, string plot_type, const int &idx = 0);
	~RTMonitor();
	bool PhotonHits(Photon &hw) const;
	void PhotonDetected(const int &w_idx, Photon *hw);
	void AddTotalPhoton(const int Nhw);
	void SetOutput(double *w_freq = NULL, int w_size = 0, string msh_file_name = "");
	void Close_oFiles(void);
	void set_file_output_off(void);
	bool Print_To_File(void) const;
	void save_to_file(void);

	dvector* Get_Output(void) const;
	int GetPhotonCount(void) const;
	dvector* get_x_data(int &x_size_out) const;
	dvector* get_y_data(int &y_size_out) const;
	double* y_data_at(const int &idx) const;
	void set_monitor(const double &xtheta, const double &xphi);
	void clean_y_data(void);

};

class predef_output {
private:
	double theta, phi;
	int stdType;
	int NumMonitors;
	string FileName;
	string *monitor_name;
	FILE *fID;
	int get_data(dvector **xdata, int * x_size, dvector **ydata, int * y_size);
	void initialize(void);

public:
	RTMonitor **std_monitor; // temporary
	predef_output(void);
	predef_output(const int &out_type,const string &xFileName);
	~predef_output(void);
	void set_object(const int &out_type, const string &xFileName);
	int link_monitor(RTMonitor *xMonitor);
	void unlink_monitors(void);
	int set_output_file(void);
	void save_to_file(void);
};

class Esource {
private:
	bool next_step;
	int iw, iPhi, iThe;
	int Nphotons, Nphi, Ntheta, Nw;
	double *phi, *theta, *w;
	doublestr Epower;
	string polarization;
	Point3D Beam_center;

	void Initiallize(void);
	void CleanUp(void);
	Point3D GeneratePosition(double xtheta,double xphi) const;
	Vector3D GeneratePolarization(void) const;
	Vector3D GenerateMomentum(double xtheta, double xphi) const;
	double *SetRange(string func, int &Nx, double C = 1);
	double *SetAngle(string str, int &Nx, string &Unit);
	void ConvertFrequency();

public:
	string wUnit, thetaUnit, phiUnit;
	Esource(int Nhw,
		string wfunc, string phifun, string thetafun, string Edistfun);
	~Esource();
	void setFrequency(string func);
	void setOrientation(string strPhi, string strTheta);
	void setEfield(string strfunc);

	int GetTotalPhotons() const;

	// Frequency data out
	double GetFrequency_val(bool printout = false) const;
	double* GetFrequency_array(bool printout = false) const;
	int GetFrequency_N(void) const;
	int GetFrequency_idx(void) const;

	// Zenith (theta) incidence angle data out
	double GetZenith_val(bool printout = false) const;
	double* GetZenith_array(void) const;	// get whole data
	int GetZenith_N(void) const;			// get array size
	int GetZenith_idx(void) const;			// get angle index

	// Azimuth (phi) incidence angle data out
	double GetAzimuth_val(bool printout = false) const;
	double* GetAzimuth_array(void) const;	// get whole data
	int GetAzimuth_N(void) const;			// get array size
	int GetAzimuth_idx(void) const;			// get angle index

	bool RunSimulation();
	bool FirstRun() const;
	Photon* MakePhoton(RTRegion **Regions, int NumRegions) const;
	bool End_spectrum(void) const;
	bool Start_spectrum(void) const;
};

class ModelBuild {
public:
	double meantest;

private:
	int nThread;
	RTMonitor **Monitor;
	RTRegion **Region;
	RTSurface **Surface;
	predef_output **std_output;
	Esource *xEsource;
	int NumSurfaces;
	int NumRegions;
	int NumMonitors;
	int NumStdOut;

	vector<InstrSet> ExtractData(ifstream &File,
		InstrSet *objinst, const string &begobj, const string &endobj, int Nmax = 1E3);
	void BuildRegions(vector<InstrSet> vRegionSet, int Ninst);
	void BuildSurfaces(vector<InstrSet> vSurfaceSet, int Ninst);
	void BuildEsource(vector<InstrSet> vEsourceSet, int Ninst);
	void BuildStdoutput(vector<InstrSet> &vMonitorSet, vector<InstrSet> vOutputSet, int Ninst);
	void BuildMonitors(vector<InstrSet> vMonitorSet, int Ninst);
	void Initiallize(void);
	void CleanUp(void);
	bool Check_SetupFile(string &setup_FileName, ifstream &File);

	void MonitorCheck(Photon *hw) const;
	void SurfaceCheck(Photon *hw) const;
	bool IsPhotonTrapped(Photon *hw, const int IsReflected) const;
	string Value(const string &arg,const vector<InstrSet> &vOjectSet,const int &iobj,const int &Narg) const;
	bool BuildObjects(const string &objName,const vector<InstrSet> &vOjectSet, InstrSet *Minobjinst,const int &Narg) const;

public :
	ModelBuild(string setupfile);
	~ModelBuild();
	void InputSetup(string setupfile);
	bool Save(void);
	bool Reset_wProperties(void);


	bool RunSimulation(void);
	int TotalPhotons(void) const;
	Photon* MakePhoton(void);
	void DetectPhotons(Photon *hw);
	void InterfaceCheck(Photon *hw) const;
	string NewPhotonState(Photon *hw) const;
	bool FirstRun() const;
	void set_wProperties(void);
};

// Other non-class functions
string CommentOut(string &str); // remove comments from text files
double RandomGen(char Type, long Seed, long *Status);
void ErrorMsg(string strError);
void GetArguments(string strArg, objArgsearch *strSearch, const int &nSearch, string strWhere = "");
string check_command_line(int argc, char* argv[]);
#endif
