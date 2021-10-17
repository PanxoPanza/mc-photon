#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <math.h>
#include "mcphoton_lib.h"

#define UNIT_FORMAT 1

#define HITSBOUNDARY	0
#define HITSSURFACE		1
#define HITSPARTICLE	2
#define ISABSORBED		3
#define HITSMONITOR		4
#define TOTALEVENTS		5


ModelBuild::ModelBuild(string setupfile) {
	// check if file exist
	Initiallize();
	InputSetup(setupfile);
}

ModelBuild::~ModelBuild() {
	CleanUp();
}

void ModelBuild::Initiallize(void) {
	nThread = 0;
	Monitor = NULL;
	Region = NULL;
	Surface = NULL;
	std_output = NULL;
	
	NumSurfaces = 0;
	NumRegions = 0;
	NumMonitors = 0;
	meantest = 0;
}

void ModelBuild::CleanUp(void) {

	if (Region != NULL) {
		for (int i = 0; i < NumRegions; i++){
			if (Region[i] != NULL) delete Region[i];
			Region[i] = NULL;
		}
		delete[] Region;
	}

	if (Surface != NULL) {
		for (int i = 0; i < NumSurfaces; i++) {
			if (Surface[i] != NULL) delete Surface[i];
			Surface[i] = NULL;
		}
		delete[] Surface;
	}

	if (std_output != NULL) {
		for (int i = 0; i < NumStdOut; i++) {
			if (std_output[i] != NULL) delete std_output[i];
			std_output[i] = NULL;
		}
		delete[] std_output;
	}

	if (Monitor != NULL) {
		for (int i = 0; i < NumMonitors; i++) {
			if (Monitor[i] != NULL) delete Monitor[i];
			Monitor[i] = NULL;
		}
		delete[] Monitor;
	}

	Initiallize();
}

bool ModelBuild::Check_SetupFile(string &setup_FileName, ifstream &File) {
	bool ask = true;
	int i_ask = 0, N_ask;

	switch (OS_SYSTEM) {
	case OS_UNI :
		N_ask = 0;
		break;
	case OS_WIN:
		N_ask = 3;
		break;
	}

	while (ask) {
		i_ask++;

		File.open(setup_FileName + ".setup");
		if (File.is_open()) {
			Log("Setup file found: %s.setup", setup_FileName.c_str());
			return true;
		}
		else {
			if (setup_FileName.empty()){
				cout << "Setup has not been specified" << endl;
				if (OS_SYSTEM == OS_UNI)
					Log("Setup has not been specified");
			}
			else {
				cout << setup_FileName << ".setup not found" << endl;
				if (OS_SYSTEM == OS_WIN)
					Log("%s.setup not found", setup_FileName.c_str());
			}

			if (i_ask <= N_ask) {
				cout << "Please enter the name of the *.setup file: ";
				cin >> setup_FileName;
			}
			else {
				ErrorMsg("Unable to run without *.setup file");
				ask = false;
			}
		}
	}
	return false;
}

void ModelBuild::InputSetup(string setupfile) {
	ifstream File;
	int Narg; // number of arguments for each class type

	
	// remove .setup extension
	size_t found = setupfile.find(".setup");
	if (found != string::npos)
		setupfile.resize(found);

	// set logfile
	switch (OS_SYSTEM) {
	case OS_UNI: // if unix, log to a .log file
		SetLogFileName(setupfile);
		break;
	case OS_WIN: // if windows log to cmd
		SetConsoleLogging();
		break;
	}

	// check for .setup file
	Check_SetupFile(setupfile, File);

	Log("Building MC-photon system...");

	InstrSet RegionSet[] = {
		{ "LABEL"		, 0, 1	, "" },
		{ "MATERIAL"	, 0, 1	, "" },
		{ "INCLUSION"	, 0, 20	, "" },
		{ "", 0, 0, "" }
	};

	InstrSet SurfaceSet[] = {
		{"LABEL"		, 0, 1	, ""},
		{"REGIONS"		, 0, 1	, "" },
		{"GEOMETRY"		, 0, 1	, "" },
		{"COATING"		, 0, 1	, "" },
		{"", 0, 0, ""}
	};

	InstrSet EsourceSet[] = {
		{ "N_PHOTONS"	, 0, 1, "" },
		{ "FREQ_RANGE"	, 0, 1, "" },
		{ "THETA"		, 0, 1, "" },
		{ "PHI"			, 0, 1, "" },
		{ "EFIELD_DIST"	, 0, 1, "" },
		{ "", 0, 0, "" }
	};

	InstrSet MonitorSet[] = {
		{"GEOMETRY"		, 0, 1, "" },
		{"FILENAME"		, 0, 1, "" },
		{"PLOT_TYPE"	, 0, 1, ""},
		{"TARGET_POINTS", 0, 1, ""},
		{"", 0, 0, ""}
	};

	InstrSet OutputSet[] = {
		{"RADIATION_PROP"	, 0, 1, "" },
		{"THREADS"			, 0, 1, "" },
		{"FILE_SUFIX"		, 0, 1, "" },
		{"", 0, 0, ""}
	};

	vector<InstrSet> vRegionSet;
	vector<InstrSet> vSurfaceSet;
	vector<InstrSet> vEsourceSet;
	vector<InstrSet> vMonitorSet;
	vector<InstrSet> vOutputSet;

	vector<string> token;
	Log("... Getting instructions from *.setup file");
	if (File.is_open()) {
		// extract data from file
		vRegionSet  = ExtractData(File, RegionSet , "REGION" , "ENDREGION"   );
		vSurfaceSet = ExtractData(File, SurfaceSet, "SURFACE", "ENDSURFACE"  );
		vMonitorSet = ExtractData(File, MonitorSet, "MONITOR", "ENDMONITOR"  );
		vEsourceSet = ExtractData(File, EsourceSet, "ESOURCE", "ENDESOURCE",1);
		vOutputSet = ExtractData(File, OutputSet  , "SETUP"  , "ENDSETUP"  ,1);

		File.close();
	}
	else {
		ErrorMsg("File: " + setupfile + ".setup not found");
	}

	/******************* Build the set of objects to run the code **************************/
	// Pass information to build regions
	Narg = (sizeof(RegionSet) / sizeof(RegionSet[0]) - 1);
	BuildRegions(vRegionSet, Narg);

	// Pass information to build surfaces
	Narg = (sizeof(SurfaceSet) / sizeof(SurfaceSet[0]) - 1);
	BuildSurfaces(vSurfaceSet, Narg);

	// Pass information to build Esource
	Narg = (sizeof(EsourceSet) / sizeof(EsourceSet[0]) - 1);
	BuildEsource(vEsourceSet, Narg);

	// Pass information to build stdoutput
	Narg = (sizeof(OutputSet) / sizeof(OutputSet[0]) - 1);
	BuildStdoutput(vMonitorSet, vOutputSet, Narg);

	// Pass information to build monitor
	Narg = (sizeof(MonitorSet) / sizeof(MonitorSet[0]) - 1);
	BuildMonitors(vMonitorSet, Narg);

	Log("... system successfully built!");
	/*************************************************************************************/
}

vector<InstrSet> ModelBuild::ExtractData(ifstream &File,
	InstrSet *objinst, const string &begobj, const string &endobj, int Nmax) {
	/* Extract the values for an object from the setup file base on a set of instructions
		- File: text file where data is stored
	//	- objinst: structure with set of names to be extracted
		- begobj: string that sets the begin of the object in text file
		- endobj: string that sets the end of the object in text file
	*/
	int objCount = 0; // count number of objects
	bool objstate = false;
	int begstr, lenstr;
	string line;
	vector<string> token;
	vector<InstrSet> objarray;

	InstrSet *objlocal;
	while (getline(File, line)) {

		CommentOut(line); // remove comments

		token = Tokenize(line); // tokenize considering " " - "\t" - "," characters

		// End of object declaration
		if (!token.empty()){ // if line is not empty
			if (objstate && !token.front().compare(endobj)) {
				for (objlocal = objinst; !objlocal->Name.empty(); objlocal++) {
					objarray.push_back(*objlocal);
					objlocal->nArgs = 0;
					objlocal->Value = "";
				}
				objstate = false; // End of object declaration
				objCount++; // count number of objects
			}

			// Object found (extract arguments)
			if (objstate) {
				for (objlocal = objinst; !objlocal->Name.empty(); objlocal++) {
					if (!token.front().compare(objlocal->Name)) { // extract line following the argument
						objlocal->nArgs++;
						if (objlocal->nArgs > objlocal->nArgsMax)
							ErrorMsg("Number of " + objlocal->Name + " declarations in " + begobj + " exceeds " + to_string(objlocal->nArgsMax));
						begstr = (int)line.find(token.at(1));
						lenstr = (int)line.size() - begstr;
						if (objlocal->nArgsMax > 1)
							objlocal->Value += "{" + line.substr(begstr, lenstr) + "} ";
						else
							objlocal->Value = line.substr(begstr, lenstr);
					}
				}
			}
			if (!token.front().compare(begobj)) objstate = true; // Start of object declaration
		}
	}

	File.clear();
	File.seekg(0, ios::beg);

	// check if the number of declarations exceeds maximum allowed
	if (objCount > Nmax) ErrorMsg("Number of " + begobj + " declarations exceeds " + to_string(Nmax));

	return objarray;
}

string ModelBuild::Value(const string &arg,const vector<InstrSet> &vOjectSet,const int &iobj,const int &Narg) const{
	for (int i = 0; i < Narg; i++){
		if (!vOjectSet.at(i + iobj*Narg).Name.compare(arg))
			return vOjectSet.at(i + iobj*Narg).Value;
	}
	return 0;
}

// Check that all the minimum instructions are defined to authorize the creation of objects
bool ModelBuild::BuildObjects(const string &objName,const vector<InstrSet> &vOjectSet,InstrSet *Minobjinst,const int &Narg) const
{
	// Check if objects have been defined
	if (vOjectSet.empty()) return false;

	// If defined, check that minimum arguments required are defined
	InstrSet *objlocal;
	int Nobj = (int)vOjectSet.size() / Narg;
	for (int i = 0; i < Nobj; i++) {
		for (objlocal = Minobjinst; !objlocal->Name.empty(); objlocal++) {
			if (Value(objlocal->Name, vOjectSet, i, Narg).empty()){
				ErrorMsg("Undefined argument " + objlocal->Name + " in " + objName);
				return false;
			}
		}
	}

	return true;
}

void ModelBuild::BuildRegions(vector<InstrSet> vRegionSet, int Ninst) {
	vector<string> token;
	string RegLabel, RegMaterial, RIncfun1, RIncfun2;
	bool setExterior = true;

	//minimum set of required arguments for each object
	InstrSet InstMin[] = {
	{ "LABEL", 0, 1, "" },
	{ "MATERIAL", 0, 1, "" },
	{ "", 0, 0, "" }
	};

	// First check if EXTERIOR has been defined
	for (int i = 0; i < vRegionSet.size() / Ninst; i++) {
		if (!Value("LABEL", vRegionSet, i, Ninst).compare("EXTERIOR")) setExterior = false;
	}
	if (setExterior) { // There si no EXTERIOR region defined
		vRegionSet.push_back(InstrSet{ "LABEL", 1, 1, "EXTERIOR" });
		vRegionSet.push_back(InstrSet{ "MATERIAL", 1, 1, "VACUUM" });
		vRegionSet.push_back(InstrSet{ "INCLUSION", 1, 1, "" });
	}

	if (!BuildObjects("REGION", vRegionSet, InstMin, Ninst))
		ErrorMsg("Undefined REGION declarations");

	// Allocate regions
	NumRegions = (int)vRegionSet.size() / Ninst;
	Region = new RTRegion*[NumRegions];

	for (int i = 0; i < NumRegions; i++) {
		// Set region's label
		token = Tokenize(Value("LABEL", vRegionSet, i, Ninst)); RegLabel = token.front();

		// Set region's material
		token = Tokenize(Value("MATERIAL", vRegionSet, i, Ninst));  RegMaterial = token.front();

		// Construct region
		Region[i] = new RTRegion(RegLabel, RegMaterial);

		// if region have inclusions
		if (!Value("INCLUSION", vRegionSet, i, Ninst).empty()) {
			Region[i]->Set_Inclusions(Value("INCLUSION", vRegionSet, i, Ninst));
		}
	}
	Log("... REGION: %i objects found", NumRegions);
}

void ModelBuild::BuildSurfaces(vector<InstrSet> vSurfaceSet, int Ninst) {
	string SRegUp, SRegLo, SGeo, SLabel;
	RTRegion *RegionUp = NULL, *RegionLo = NULL;
	vector<string> token;
	bool FoundRegUp, FoundRegLo;

	// ************  Check instructions before constructing object ****************
	//		Minimum set of required arguments per object
	InstrSet InstMin[] = {
	{ "LABEL", 0, 1, "" },
	{ "REGIONS", 0, 1, "" },
	{ "GEOMETRY", 0, 1, "" },
	{ "", 0, 0, "" }
	};

	if (!BuildObjects("SURFACE", vSurfaceSet, InstMin, Ninst)) // If false gives error
		ErrorMsg("Undefined SURFACE declarations");
	//*********************************************************************************

	// Allocate surfaces
	NumSurfaces = (int)vSurfaceSet.size() / Ninst;
	Surface = new RTSurface*[NumSurfaces];

	for (int i = 0; i < NumSurfaces; i++) {
		FoundRegUp = false; FoundRegLo = false;

		// Get surface label
		token = Tokenize(Value("LABEL", vSurfaceSet, i, Ninst));
		SLabel = token.front();

		// Get upper and lower regions
		token = Tokenize(Value("REGIONS", vSurfaceSet, i, Ninst));
		SRegUp = token.at(0);
		SRegLo = token.at(1);
		for (int j = 0; j < NumRegions; j++) {
			if (!Region[j]->Label.compare(SRegUp)) {
				RegionUp = Region[j];
				FoundRegUp = true;
			}
			if (!Region[j]->Label.compare(SRegLo)) {
				RegionLo = Region[j];
				FoundRegLo = true;
			}
		}

		if (!(FoundRegUp && FoundRegLo))
			ErrorMsg("Undefined regions for Surface " + SLabel);

		// Get the geometry function
		SGeo = Value("GEOMETRY", vSurfaceSet, i, Ninst);
		Surface[i] = new RTSurface(RegionUp, RegionLo, SGeo, SLabel, i);

		// If surface has multilayer coatings
		if (!Value("COATING", vSurfaceSet, i, Ninst).empty()) {
			Surface[i]->SetCoating(Value("COATING", vSurfaceSet, i, Ninst));
		}
	}
	Log("... SURFACE: %i objects found", NumSurfaces);
}

void ModelBuild::BuildEsource(vector<InstrSet> vEsourceSet, int Ninst) {
	int ii;
	vector<string> token;
	string RegLabel, RegMaterial, RIncfun1, RIncfun2;
	bool setExterior = true;

	// ************  Check instructions before constructing object ****************
	//		Minimum set of required arguments per object
	InstrSet InstMin[] = {
		{ "N_PHOTONS", 0, 1, "" },
		{ "FREQ_RANGE", 0, 1, "" },
		{ "THETA", 0, 1, "" },
		{ "PHI", 0, 1, "" },
		{ "EFIELD_DIST", 0, 1, "" },
		{ "", 0, 1, "" }
	};
	//		 if false gives error
	if (!BuildObjects("ESOURCE", vEsourceSet, InstMin, Ninst))
		ErrorMsg("Undefined ESOURCE declarations");
	//*********************************************************************************

	ii = 0;
	// Get number of photon
	token = Tokenize(Value("N_PHOTONS", vEsourceSet, 0, Ninst));
	int Nphotons = (int)stod(token.front());
	Log("... SOURCE: %i photons per frequency", Nphotons);

	// Get frequency-function
	string wfun = Value("FREQ_RANGE", vEsourceSet, 0, Ninst);
	Log("... SOURCE: Frequency range, '%s'", wfun.c_str());

	// Get theta function
	string thetafun = Value("THETA", vEsourceSet, 0, Ninst);
	Log("... SOURCE: Theta angle range, '%s'", thetafun.c_str());

	// Get theta function
	string phifun = Value("PHI", vEsourceSet, 0, Ninst);
	Log("... SOURCE: Phi angle range, '%s'", thetafun.c_str());

	// Get theta function
	string Edistfun = Value("EFIELD_DIST", vEsourceSet, 0, Ninst);
	Log("... SOURCE: Field distribution, '%s'", Edistfun.c_str());

	// Construct region
	xEsource = new Esource(Nphotons, wfun, phifun, thetafun, Edistfun);

	Log("... SOURCE: object successfully built");
}

void ModelBuild::BuildStdoutput(vector<InstrSet> &vMonitorSet, vector<InstrSet> vOutputSet, int Ninst) {
	int ii;
	NumStdOut = 0;
	vector<string> out_file;
	vector<int> out_type;
	vector<string> token;
	InstrSet InstMin[] = {
		{"", 0, 0, ""}
	};

	// ************  Check instructions before constructing object ****************
	if (!BuildObjects("SETUP", vOutputSet, InstMin, Ninst)) return;
	//*********************************************************************************

	ii = 0;
	if (!Value("RADIATION_PROP", vOutputSet, 0, Ninst).empty()) { // power radiation output
		token = Tokenize(Value("RADIATION_PROP", vOutputSet, 0, Ninst));
		out_file.push_back(token.at(0));
		out_type.push_back(RAD_PROPERTIES);

		InstrSet MonitorSet[] = {
		{"GEOMETRY", 0, 1, "" },
		{"FILENAME", 0, 1, "" },
		{"REGIONS" , 0, 1, "" },
		{"PLOT_TYPE", 0, 1, ""},
		{"TARGET_POINTS", 0, 1, ""},
		};

		// set upper monitor for reflectivity
		MonitorSet[0] = { "GEOMETRY", 1, 1, "BOX_Z-OPEN(1E8,1E8,1E5) CENTER(0,0,+0.5E5)" };
		MonitorSet[1] = { "FILENAME", 1, 1, "R_" + out_file.at(NumStdOut) };
		MonitorSet[2] = { "REGIONS", 1,  1, "EXTERIOR EXTERIOR" };
		MonitorSet[3] = { "PLOT_TYPE", 1,  1, "RADIATION PROPERTIES" };
		MonitorSet[4] = { "TARGET_POINTS", 1, 1,  "" };
		vMonitorSet.push_back(MonitorSet[0]);
		vMonitorSet.push_back(MonitorSet[1]);
		vMonitorSet.push_back(MonitorSet[2]);
		vMonitorSet.push_back(MonitorSet[3]);
		vMonitorSet.push_back(MonitorSet[4]);

		// set lower monitor for transmissivity
		MonitorSet[0] = { "GEOMETRY", 1, 1, "BOX_Z+OPEN(1E8,1E8,1E5) CENTER(0,0,-0.5E5)" };
		MonitorSet[1] = { "FILENAME", 1, 1, "T_" + out_file.at(NumStdOut) };
		MonitorSet[2] = { "REGIONS", 1,  1, "EXTERIOR EXTERIOR" };
		MonitorSet[3] = { "PLOT_TYPE", 1,  1, "RADIATION PROPERTIES" };
		MonitorSet[4] = { "TARGET_POINTS", 1, 1, "" };
		vMonitorSet.push_back(MonitorSet[0]);
		vMonitorSet.push_back(MonitorSet[1]);
		vMonitorSet.push_back(MonitorSet[2]);
		vMonitorSet.push_back(MonitorSet[3]);
		vMonitorSet.push_back(MonitorSet[4]);

		NumStdOut++;;
	}

	// Check for number of threads declared on the script
	if (!Value("THREADS", vOutputSet, 0, Ninst).empty()){ // if the number is fixed
		nThread = (int)stod(Value("THREADS", vOutputSet, 0, Ninst));
		if (nThread > omp_get_max_threads()){
			Log("... SETUP: number of threads exceeds maximum available");
			nThread = omp_get_max_threads();
			}
		}
	else	// no declaration, use maximum number of threads
		nThread = omp_get_max_threads();
	
	omp_set_num_threads(nThread);
	Log("... SETUP: setting OpenMP multithreading (%i threads)", nThread);
	
	file_sufix = Value("FILE_SUFIX", vOutputSet, 0, Ninst);

	Log("... SETUP: %i STANDARD OUTPUT declarations found", NumStdOut);
	if (NumStdOut > 0) {
		std_output = new predef_output*[NumStdOut];
		for (int i = 0; i < NumStdOut; i++)
			std_output[i] = new predef_output (out_type.at(i), out_file.at(i));
	}

	Log("... SETUP: object successfully built");
}

void ModelBuild::BuildMonitors(vector<InstrSet> vMonitorSet, int Ninst) {
	string MonGeo, oFile, plot_type, xy_target_file;
	vector<string> token;
	string MonRegUp, MonRegLo;
	bool FoundRegUp, FoundRegLo;
	int TotalPhotons = xEsource->GetTotalPhotons();
	RTRegion *RegionUp = NULL, *RegionLo = NULL;

	// Extract frequency from the source
	double *w_freq = xEsource->GetWavelength_array();
	int w_size = xEsource->GetWavelength_N();

	// ************  Check instructions before constructing object ****************
	//		minium set of required arguments per object
	InstrSet InstMin[] = {
		{"GEOMETRY", 0, 1, "" },
		{"FILENAME", 0, 1, "" },
		{"PLOT_TYPE", 0, 1, ""},
		{"", 0, 0, ""}
	};
	//		 if false gives error
	if (!BuildObjects("MONITOR", vMonitorSet, InstMin, Ninst)) ErrorMsg("Undefined MONITOR declarations");
	//*********************************************************************************

	// Allocate monitors
	NumMonitors = (int)vMonitorSet.size() / Ninst;
	Monitor = new RTMonitor*[NumMonitors];

	for (int i = 0; i < NumMonitors; i++) {

		// Set geometry
		MonGeo = Value("GEOMETRY", vMonitorSet, i, Ninst);

		// set output file
		token = Tokenize(Value("FILENAME", vMonitorSet, i, Ninst)); 
		oFile = token.front() + '_' + file_sufix;
		
		// get plot type
		plot_type = Value("PLOT_TYPE", vMonitorSet, i, Ninst);

		// get xy_target points file
		xy_target_file = Value("TARGET_POINTS", vMonitorSet, i, Ninst);

		// construct monitor
		Monitor[i] = new RTMonitor(MonGeo, oFile, TotalPhotons, plot_type, i);

		// check if monitor is linked to a standard output object
		for (int j = 0; j < NumStdOut; j++)
			if (std_output[j]->link_monitor(Monitor[i])) break;

		Monitor[i]->SetOutput(w_freq, w_size, xy_target_file, nThread);
	}

	Log("... MONITOR: %i objects found", NumMonitors - 2 * NumStdOut);
}

void ModelBuild::InterfaceCheck(Photon *hw) const{

	SurfaceCheck(hw);// check if photon hits surfaces and gets closest distance
	MonitorCheck(hw);// check if photon hits monitors and gets closest distance
}

void ModelBuild::SurfaceCheck(Photon *hw) const{
	for (int i = 0; i < NumSurfaces; i++)
		Surface[i]->PhotonHits(*hw);
}

void ModelBuild::MonitorCheck(Photon *hw) const{
	for (int i = 0; i < NumMonitors; i++)
		Monitor[i]->PhotonHits(*hw);
}

string ModelBuild::NewPhotonState(Photon *hw) const {
	//int NumStates = TOTALEVENTS;
	string State = "NULL";
	int FresnelState = 0;
	//int SurfaceID;
	//int MonitorID;

	double *Length = new double[TOTALEVENTS];
	Length[HITSBOUNDARY] = hw->GetBoundaryDistace();	// Distance to exit the domain
	Length[HITSSURFACE] = hw->GetSurfaceDistance();	// Distance to hit a surface
	Length[ISABSORBED] = hw->GetExtLength();			// distance to be absorbed by the material
	Length[HITSPARTICLE] = hw->GetParticleLength();	// Distance to hit a particle
	Length[HITSMONITOR] = hw->GetMonitorDistance();	// Distance to hit a particle

	switch (indexMin(Length, TOTALEVENTS)) { // Get index of minimum length
	case HITSBOUNDARY: // Photon Hits Boundary
		State = "Photon Hits Boundary";
		hw->PhotonKill();
		break;

	case HITSSURFACE: // Photon Hits Interface
		State = "Photon Hits Surface: " + Surface[hw->GetSurfaceID()]->Label;
		FresnelState = Surface[hw->GetSurfaceID()]->RefractPhoton(*hw);
		if (FresnelState == 1) {
			State += " (Reflected)";
		}
		else if (FresnelState == -1) {
			State += " (Transmitted)";
		}
		else if (FresnelState == 0) {
			State += " (Absorbed)";
		}
		break;

	case HITSPARTICLE: // Photon hits a particle
		State = "Photon hits a particle " + hw->GetParticleLabel();
		if (hw->HitsParticle()) {
			State += "(Scattered)";
			FresnelState = 1;
		}
		else
			State += "(Absorbed)";
		break;

	case ISABSORBED: // Photon is absorbed by region
		State = "Photon absorbed in: " + hw->GetRegion()->Label;
		hw->PhotonKill();
		break;

	case HITSMONITOR: // Photon hits a monitor
		int iw = xEsource->GetWavelength_idx();

		State = "Photon hits Monitor: " + Monitor[hw->GetMonitorID()]->Label;
		Monitor[hw->GetMonitorID()]->PhotonDetected(iw, hw);
		break;
	}
	delete[] Length;

	// Kill photon if goes internal reflection
	if (IsPhotonTrapped(hw, FresnelState)) {
		hw->PhotonKill();
		State = State + ", Multiple reflections";
	}

	return State;
}

void ModelBuild::DetectPhotons(Photon *hw){
	int NumStates = 4;
	int iw = xEsource->GetWavelength_idx();
	string State = "NULL";

	// Determine if photon is closer to monitor[i]
	double *Length = new double[NumStates];
	Length[0] = hw->GetMonitorDistance();	// Distance to hit a monitor
	Length[1] = hw->GetSurfaceDistance();	// Distance to hit a surface
	Length[2] = hw->GetExtLength();			// distance to be absorbed by the material
	Length[3] = hw->GetParticleLength();	// Distance to hit a particle

	// Get index of minimum length
	if (indexMin(Length, NumStates) == 0) {
		int i = hw->GetMonitorID();
		Monitor[i]->PhotonDetected(iw,hw);
	}
	delete[] Length;
}

bool ModelBuild::IsPhotonTrapped(Photon *hw, const int IsReflected) const
/*****************************************************************************************
	Detects whether the photon is trapped inside a material due
	to multiple internal reflections
*****************************************************************************************/
{
	if (IsReflected == 1) { // If photon reflected
		string newRegionLabel = hw->GetRegion()->Label;
		if (hw->NumReflects > 0 && newRegionLabel.compare(hw->oldRegionLabel) == 0) {
			hw->NumReflects++;
			hw->oldRegionLabel = newRegionLabel; // record region
			if (hw->NumReflects > MAX_REFLECTION) { // consecutive internal reflections or scattering events
				hw->NumReflects = 0;
				return true;
			}

		}
		else {	// Photon goes first reflection
			hw->NumReflects = 1;
			hw->oldRegionLabel = newRegionLabel; // record region
			return false;
		}
	}
	else {			// Photon is not reflected
		hw->NumReflects = 0;
		hw->oldRegionLabel = hw->GetRegion()->Label;
	}
	return false;
}

Photon* ModelBuild::MakePhoton(void)  {
	//meantest -= log(RandomNum);
	//meantest++;
	return xEsource->MakePhoton(Region,NumRegions);
}

bool ModelBuild::RunSimulation(void) {
	bool run_simulation = xEsource->RunSimulation();

	if (!run_simulation)
		Log("done with calculations ...");

	return run_simulation;
}

int ModelBuild::TotalPhotons(void) const {
	return xEsource->GetTotalPhotons();
}

bool ModelBuild::FirstRun(void) const {
	return xEsource->FirstRun();
}

void ModelBuild::set_wProperties(void) {
	// Erase all wProperties
	if (!Reset_wProperties()) ErrorMsg("Unable to reset wProperties");

	double xtheta = xEsource->GetZenith_val(UNIT_FORMAT);
	double xphi = xEsource->GetAzimuth_val(UNIT_FORMAT);
	double w = xEsource->GetWavelength_val();
	double lambda = xEsource->GetWavelength_val(UNIT_FORMAT);
	const char *Unitlambda = xEsource->wUnit.c_str();

	// Set angle of incidence
	if (xEsource->Start_spectrum()) {
		Log("Setting angle of incidence (theta phi) = (%6.2f, %6.2f) deg", xtheta, xphi);
		for (int i = 0; i < NumMonitors; i++)
			Monitor[i]->set_monitor(xtheta, xphi);
	}

	// Set frequency values to regions and surfaces
	for (int i = 0; i < NumRegions; i++)
		if (!Region[i]->setProperties(w))
			ErrorMsg("Error: Failed to set properties in Region: "
				+ Region[i]->Label + " at " + to_string(lambda));

	for (int i = 0; i < NumSurfaces; i++)
		if (!Surface[i]->set_wLayers(w))
			ErrorMsg("Error: Failed to set thinfilm properties in Surface: "
				+ Surface[i]->Label + " at " + to_string(lambda));

	// output message
	Log("... Computing for lambda = %7.3f %-s", lambda, Unitlambda);

	/*cout << "lambda = " + to_string(lambda)
	+ " eps = "
	+ to_string(real(Region[0]->get_wEps())) + " + "
	+ to_string(imag(Region[0]->get_wEps())) + "i" << endl;
	*/
	meantest = 0.0;
}

// Reset saved properties at frequency w (see "set_wProperties")
bool ModelBuild::Reset_wProperties(void) {

	// Reset Regions
	for (int i = 0; i < NumRegions; i++)
		Region[i]->reset_wProperties();

	// Reset Surfaces
	for (int i = 0; i < NumSurfaces; i++)
		Surface[i]->reset_wProperties();

	return true;
}


bool ModelBuild::Save(void) {

	if (!xEsource->End_spectrum())
		return false;
	
	Log("... Writting results to file ...");
	// write results to file
	for (int i = 0; i < NumMonitors; i++)
		if (Monitor[i]->Print_To_File())
			Monitor[i]->save_to_file();

	for (int i = 0; i < NumStdOut; i++)
		std_output[i]->save_to_file();

	// clean monitors
	for (int i = 0; i < NumMonitors; i++)
		Monitor[i]->clean_y_data();

	return true;
}
