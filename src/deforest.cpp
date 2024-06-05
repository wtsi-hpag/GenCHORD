#include "../libs/JSL/JSL.h"
// #define GNUPLOT_NO_TIDY
#include "Utility/basicFunctions.h"
#include "DataFrame/data.h"
#include "ParameterInference/GlobalInference.h"
#include "Utility/plotting.h"
// #include "Utility/logFactorial.h"

#include "HarmonicTree/GetHarmonics.h"
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

	Log("Preparing probability models:\n")
	int Kmax = 2*d.Mean + d.Deviation;
	int qmax = settings.Qmax;
	auto model = Models::NegativeBinomial();

	NormaliseModel(model,d,settings, Kmax,qmax);

	// model.SetSignalParameters(11.60,5.33);
	// model.SetNoiseParameters(10.1452,4.7,0.30);
	
	Kmax = d.maxK;
	model.SetDimensionality(Kmax,qmax);
	model.SetGrids();


	Log("Beginning network navigation" << std::endl;)
	auto path = GetHarmonics(d,settings,model);
	JSL::gnuplot gp;

	settings.MemorySmoothing = 0.999;
	// Data d2(settings);

	basicPlot(gp,d,0);
	TransitionPlot(gp,d,path,"test");
	gp.SetPersistence(true);
	gp.Show();

	Log("Deforest routine completed. Have a nice day.\n\n")
}