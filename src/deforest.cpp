// #define GNUPLOT_NO_TIDY
#include "../libs/JSL/JSL.h"
// #include "Utility/plotting.h"
#include "Utility/basicFunctions.h"
#include "DataFrame/data.h"
#include "Utility/logFactorial.h"

// #include "HarmonicTree/GetHarmonics.h"
#include "settings.h"
#include <chrono>
using namespace std::chrono;
LogFactorial LogFac; //initialise the functor object
int nPlot = 8;

bool globalRegress = true;
int main(int argc, char**argv)
{
	//read the command line arguments first
	JSL::Argument<std::string> configFile("__none__","config",argc,argv);
	Settings settings;
	if (configFile.Value == "__none__")
	{
		settings.Initialise(argc,argv);
	}
	else
	{
		settings.Initialise(configFile.Value);
	}

	//some pretty introductory text
	globalVerbose = !settings.Quiet;
	Log("==========================================\n\tCoverage Deforesting\n==========================================" << std::endl);
	

	//load the data file -- either from file, or from a pipe
	Data d(settings);
	// if (JSL::PipedInputFound())
	// {
	// 	// d = Data(settings.DataThinning,settings.TargetChromosome,0.99);
	// }
	// else
	// {
	// 	d = Data(settings.DataFile,settings.DataThinning,settings.TargetChromosome,0.99);
	// }

	
	
	// JSL::gnuplot gp;
	
	// auto path = GetHarmonics(d,settings,gp);


	Log("Deforest routine completed. Have a nice day.\n\n")
}