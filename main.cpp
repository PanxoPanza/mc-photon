#define _USE_MATH_DEFINES
#include "MCRT_library.h"
#define RandNum0 rand() % 100 // Random generator for seeds

Photon *hw;
ofstream debug_file;
string ModelBuild_file;
double w;
bool debug = false;
/************************************* Main Code **********************************/
int main(int argc, char* argv[])
{
	string PhotonState;
	string debug_stream = "";

	ModelBuild_file = check_command_line(argc, argv);
	//ModelBuild_file = "H2O_catalyst_film.setup";
	ModelBuild Setup(ModelBuild_file);

	if (debug) debug_file.open("debug_file.txt");
	if (debug) omp_set_num_threads(1);

	srand(10);
	int NT = omp_get_num_threads();
	while (Setup.RunSimulation()) {


		// Set material properties for a given frequency
		Setup.set_wProperties();

		/*
		int count[NT], totalcount = 0;
		for(int ii; ii < NT; ii++)
				count[ii] = 0;


		//Log("parameters set, total = %i, number of processors %i", totalcount,  omp_get_num_threads());
		*/

		// Run parallel calculation for all photons
		#pragma omp parallel default(none), private(hw, PhotonState), shared(Setup, debug_stream, debug, debug_file)
		{

			// get seed for each thread on first run
			if (Setup.FirstRun()) RandomGen(0, RandNum0, NULL);

			#pragma omp for
			for (int ihw = 0; ihw < Setup.TotalPhotons(); ihw++)
			{
				int nt = omp_get_thread_num();
				hw = Setup.MakePhoton();
				Setup.meantest -= log((double) RandomNum);

				if (debug) debug_file
				<< "********************** New Photon **********************\n";

				while (hw->IsAlive()) {
					// Clean photon interface check
					hw->InitializeLength();

					// Check distance to all interfaces
					if (debug) debug_file
						<< "Interface Check: " << endl;

					Setup.InterfaceCheck(hw);
					if (debug) debug_file
						<< hw->Print();

					// Get a new photon state
					PhotonState = Setup.NewPhotonState(hw);

					if (debug) debug_file
						<< "	New Photon State: " + PhotonState + "\n";
				}
				delete hw;
			}
		}
		/*
		for (int ii; ii < NT; ii++){
			totalcount += count[ii];
			Log("count[%i] = %.1f", ii, count[ii]);
		}*/

		// Print results to a file
		Log("mean of log(rand) %.5f",Setup.meantest/Setup.TotalPhotons());
		Setup.Save();
	}
	return 0;
}
