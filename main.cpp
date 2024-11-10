#define _USE_MATH_DEFINES
#include "mcphoton_lib.h"
#include <cstdlib> // for sprintf
#define RandNum0 rand() % 100 // Random generator for seeds

ofstream debug_file;
double w;
bool debug;
/************************************* Main Code **********************************/
int main(int argc, char* argv[])
{
	//string debug_stream = "";

	string ModelBuild_file = check_command_line(argc, argv);
	ModelBuild Setup(ModelBuild_file);
	debug = Setup.debug_mode;

	if (debug) {
		debug_file.open("debug_file.txt");
		omp_set_num_threads(1);
	}

	srand(10);
	while (Setup.RunSimulation()) {


		// Set material properties for a given frequency
		Setup.set_wProperties();

		// Run parallel calculation for all photons
		#pragma omp parallel default(none), shared(Setup, debug, debug_file)
		{

			// get seed for each thread on first run
			if (Setup.FirstRun()) RandomGen(0, RandNum0, NULL);

			#pragma omp for schedule(dynamic,1)
			for (int ihw = 0; ihw < Setup.TotalPhotons(); ihw++)
			{
				Photon *hw = Setup.MakePhoton();

				if (debug) debug_file
				<< "********************** New Photon **********************\n";

				while (hw->IsAlive()) {
					// Clean photon interface check
					hw->InitializeLength();

					// Check distance to all interfaces
					if (debug) debug_file << "Interface Check: " << endl;
					Setup.InterfaceCheck(hw);
					
					if (debug) debug_file << hw->Print();

					// Get a new photon state
					string PhotonState = Setup.NewPhotonState(hw);

					if (debug) debug_file
						<< "	New Photon State: " + PhotonState + "\n";
				}
				delete hw;
			}
		}
		Setup.Save();
	}
	return 0;
}
