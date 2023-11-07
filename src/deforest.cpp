#define GNUPLOT_NO_TIDY
#include "../libs/JSL/JSL.h"
#include "Utility/plotting.h"
#include "Utility/basicFunctions.h"
#include "data.h"
// #include "basicSmooth.h"
#include "Utility/logFactorial.h"
#include "HarmonicBayes/harmonicFit.h"

#include "HarmonicTree/tree.h"
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
	
	Data d;
	if (!settings.PlotOnly)
	{
		d = Data(settings.DataFile,settings.DataThinning,settings.TargetChromosome,settings.MemorySmoothing);

		TreeAssign(d,settings);
		// HarmonicFit(d,settings);
		// settings.DataFile = settings.OutputName;
		// OutputPlot(d,settings);
	}
	else
	{
		if (settings.ComparePlot == "__none__")
		{
			OutputPlot(d,settings);
		}
		else
		{
			ComparisonPlots(settings);
		}
	}

	Log("Deforest routine completed. Have a nice day.\n\n")
}