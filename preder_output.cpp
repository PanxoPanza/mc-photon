#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include "mcphoton_lib.h"

// a class for standard analysis. It uses the data gathered by the monitors.
// warning!: this class does not own the monitors and therefore it may not delete the monitors, 
//			but rather unlink the pointers (point to NULL) 
predef_output::predef_output(void) {
	initialize();
}

predef_output::predef_output(const int &out_type, const string &xFileName) {

	initialize();
	set_object(out_type, xFileName);
}

void predef_output::set_object(const int &out_type, const string &xFileName) {
	stdType = out_type;
	FileName = xFileName;
	fID = NULL;
	switch (stdType) {
	case RAD_PROPERTIES:
		NumMonitors = 2;

		// define monitors' name
		monitor_name = new string[NumMonitors];
		monitor_name[0] = "R_" + FileName;
		monitor_name[1] = "T_" + FileName;
		break;

	default:
		ErrorMsg("SETUP, Unknown standard output type");
		break;
	}

	// create monitor pointer's array
	std_monitor = new RTMonitor *[NumMonitors];
	unlink_monitors();
	set_output_file();
}

predef_output::~predef_output(void) {
	if (std_monitor != NULL) { 
		unlink_monitors();
		delete[] std_monitor;
	}

	if (monitor_name != NULL) {
		for (int i = 0; i < NumMonitors; i++) 
			monitor_name[i].clear();
		delete[] monitor_name;
	}
	if (fID != NULL) fclose(fID);
	initialize();
}


void predef_output::initialize(void) {
	fID = NULL;
	std_monitor = NULL;
	monitor_name = NULL;
	NumMonitors = 0;
	FileName.clear();
}
void predef_output::unlink_monitors(void) {
	for (int i = 0; i < NumMonitors; i++)
		if (std_monitor[i] != NULL) std_monitor[i] = NULL;
}

// link an external monitor pointer with the internal monitors
// returns 1 if successfull
int predef_output::link_monitor(RTMonitor *xMonitor) {
	int found_monitor = 0; // state set to 0 by default

	// check if monitor array has been allocated
	if (NumMonitors == 0 || std_monitor == NULL)
		ErrorMsg("Undefined monitors in predef_output object " + FileName);
	
	for (int idx_mon = 0; idx_mon < NumMonitors; idx_mon++) {
		if (std_monitor[idx_mon] == NULL &&						// Check if monitor pointer has not being defined
			!monitor_name[idx_mon].compare(xMonitor->Label)) {  // Looks for monitor label
			
			//link monitor
			std_monitor[idx_mon] = xMonitor;
			
			// monitor's own ouput is turned off as output is set inside the class
			std_monitor[idx_mon]->set_file_output_off();
			
			// monitor found
			found_monitor = 1;
			break;
		}
	}

	return found_monitor;
}

// create output files. Returns 1 if successfull
int predef_output::set_output_file(void) {
	string oFile;

	switch (stdType) {
	case RAD_PROPERTIES:
		oFile = FileName + ".rad_prop";
		if (fID != NULL) return 0;
		fID = fopen(oFile.c_str(), "w");
		fprintf(fID, "# Radiation properties file\n");
		fprintf(fID, "# File columns:\n");
		fprintf(fID, "# \t 1 Incidence zenith angle (deg) \n");
		fprintf(fID, "# \t 2 Incidence azimuth angle (deg) \n");
		fprintf(fID, "# \t 3 Frequency (rad/s)\n");
		fprintf(fID, "# \t 4 Reflectivity\n");
		fprintf(fID, "# \t 5 Transmittance\n");
		fprintf(fID, "# \t 6 Reflected flux\n");
		fprintf(fID, "# \t 7 Transmitted flux\n");
		fprintf(fID, "# \t 8 Mean distance (um) traveled by reflected photons\n");
		fprintf(fID, "# \t 9 Mean distance (um) traveled by transmitted photons \n");
		fprintf(fID, "# \n");
		break;
	default :
		ErrorMsg("SETUP, Unknown standard output type");
		break;
	}

	return 1;
}

void predef_output::save_to_file(void) {
	dvector **x_data, **y_data;

	int *x_size, *y_size;
	if (fID == NULL) ErrorMsg("Monitor file not set");

	// Allocate data storage variables
	x_size = new int[NumMonitors];
	y_size = new int[NumMonitors];
	x_data = new dvector*[NumMonitors];
	y_data = new dvector*[NumMonitors];
	
	// extract the data from linked monitors
	//fprintf(fID, "# Date of calculation:\t");
	//fprintf(fID, ctime(&end_time)); // get actual date of calculation
	//fprintf(fID, "# \n");
	switch (stdType) {
	case RAD_PROPERTIES:
		if(!get_data(x_data, x_size, y_data, y_size))// extract data collected by monitors
			ErrorMsg("Unable to extract data in standard output object:\n" + FileName); 

		for (int iw = 0; iw < x_data[0][0].size(); iw++) {
			// print x values

			fprintf(fID, "%-6.2f\t", std_monitor[0]->theta); // Zenith
			fprintf(fID, "%-6.2f\t", std_monitor[0]->phi); // Azimuth
			fprintf(fID, "%-7.4e\t",  x_data[0][0].at(iw)); // frequency
			
			// print y values
			double Rdetect = y_data[0][0].at(iw) == 0 ? 1E-5 : y_data[0][0].at(iw);
			double Tdetect = y_data[1][0].at(iw) == 0 ? 1E-5 : y_data[1][0].at(iw);
			fprintf(fID, "%-9.6e\t", y_data[0][0].at(iw)); // reflection
			fprintf(fID, "%-9.6e\t", y_data[1][0].at(iw)); // transmission
			fprintf(fID, "%-9.6e\t", y_data[0][1].at(iw)/ Rdetect); // reflected flux
			fprintf(fID, "%-9.6e\t", y_data[1][1].at(iw)/ Tdetect); // transmited flux
			fprintf(fID, "%9.6e\t",  y_data[0][2].at(iw)/ Rdetect); // Mean path from R photons
			fprintf(fID, "%9.6e\n",  y_data[1][2].at(iw)/ Tdetect); // Mean path from T photons
		}
		break;

	default:
		ErrorMsg("SETUP, Unknown standard output type");
		break;
	}

	// Free x_data and y_data pointers and delete the arrays
	for (int i = 0; i < NumMonitors; i++) {
		x_data[i] = NULL; // unlink x_data[i] pointer without deleting it
		y_data[i] = NULL; // unlink y_data[i] pointer without deleting it
	}
	delete[] x_data;
	delete[] y_data;
}

// this function creates an array of pointers to the data inside the monitors
int predef_output::get_data(dvector **xdata, int * x_size, dvector **ydata, int * y_size) {

	for (int i = 0; i < NumMonitors; i++) {
		xdata[i] = std_monitor[i]->get_x_data(x_size[i]);
		ydata[i] = std_monitor[i]->get_y_data(y_size[i]);
	}
	return 1;
}
