// #define GNUPLOT_NO_TIDY
#include "JSL.h"
#include "plotting.h"
#include "basicFunctions.h"
#include "data.h"
// #include "basicSmooth.h"
#include "logFactorial.h"
#include "harmonicFit.h"
#include "settings.h"
#include <chrono>
using namespace std::chrono;
LogFactorial LogFac; //initialise the functor object
int nPlot = 8;

bool globalRegress = true;
int main(int argc, char**argv)
{
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
	globalVerbose = !settings.Quiet;
	Log("==========================================\n\tCoverage Deforesting\n==========================================" << std::endl);
	

	if (!settings.PlotOnly)
	{
		Data d(settings.DataFile,settings.DataThinning,settings.TargetChromosome,settings.MemorySmoothing);
		HarmonicFit(d,settings);
		settings.DataFile = settings.OutputName;
		OutputPlot(settings);
	}
	else
	{
		if (settings.ComparePlot == "__none__")
		{
			OutputPlot(settings);
		}
		else
		{
			ComparisonPlots(settings);
		}
	}

	Log("Deforest routine completed. Have a nice day.\n\n")
}