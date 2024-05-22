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

int main(int argc, char**argv)
{
	//read the command line arguments first
	JSL::Argument<std::string> configFile("__none__","config",argc,argv);
	Settings settings(configFile,argc,argv);
	
	//some pretty introductory text
	globalVerbose = !settings.Quiet;
	Log("==========================================\n\tCoverage Deforesting\n==========================================" << std::endl);
	

	//load the data file -- either from file, or from a pipe
	Data d(settings);

	
	
	JSL::gnuplot gp;

	gp.Plot(d.Chromosomes[0].Idx,d.Chromosomes[0].Counts);
	gp.Show();
	
	// auto path = GetHarmonics(d,settings,gp);


	Log("Deforest routine completed. Have a nice day.\n\n")
}