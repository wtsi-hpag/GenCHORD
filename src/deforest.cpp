#include "../libs/JSL/JSL.h"
// #define GNUPLOT_NO_TIDY
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
	Data d(settings);
	int Kmax = std::min(150,d.maxK);
	int qmax = 10;
	// int Kmax = 200;
	auto model = Models::NegativeBinomial(25,5,-3,100,10);
	std::vector<double> ws = {0.1,0.05,0.65,0.2};

	// auto d = model.Draw(100000,ws,Kmax);

	// JSL::gnuplot;
	// gp.Plot(d.Chromosomes[0])
	// JSL::gnuplot gp;
	
	// settings.DataThinning = 1e5;
	// Data d2(settings);
	// gp.Plot(d2.Chromosomes[0].Idx,d2.Chromosomes[0].Counts);
	// gp.Show();
	GlobalInference(model,d,settings, Kmax,qmax);
	
	// auto path = GetHarmonics(d,settings,gp);


	Log("Deforest routine completed. Have a nice day.\n\n")
}