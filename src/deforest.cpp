// #define GNUPLOT_NO_TIDY
#include "../libs/JSL/JSL.h"
#include "Utility/basicFunctions.h"
#include "DataFrame/data.h"
#include "ParameterInference/GlobalInference.h"
// #include "Utility/plotting.h"
// #include "Utility/logFactorial.h"

#include "settings.h"
#include <chrono>
using namespace std::chrono;
// LogFactorial LogFac; //initialise the functor object

int main(int argc, char**argv)
{
	//read the command line arguments first
	JSL::Argument<std::string> configFile("__none__","config",argc,argv);
	Settings settings(configFile,argc,argv);
	
	//some pretty introductory text
	globalVerbose = !settings.Quiet;
	Log("==========================================\n\tCoverage Deforesting\n==========================================" << std::endl);
	

	//load the data file -- either from file, or from a pipe
	// Data d(settings);

	auto model = Models::Gaussian(25,2,-1,50,160);
	
	int Kmax = 150;
	std::vector<double> ws = {0.1,0.6,0.1,0.2};
	model.Normalise(Kmax,ws);

	auto d = model.Draw(2000000,ws,Kmax);


	GlobalInference(model,d,settings,Kmax,14);
	// JSL::gnuplot gp;

	// gp.Plot(d.Chromosomes[0].Idx,d.Chromosomes[0].Counts);
	// gp.Show();
	
	// auto path = GetHarmonics(d,settings,gp);


	Log("Deforest routine completed. Have a nice day.\n\n")
}