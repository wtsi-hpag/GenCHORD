#define GNUPLOT_NO_TIDY
#include "../libs/JSL/JSL.h"
#include "Utility/plotting.h"
#include "Utility/basicFunctions.h"
#include "data.h"
#include "Utility/logFactorial.h"

#include "HarmonicTree/GetHarmonics.h"
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
	
	JSL::gnuplot gp;
	d = Data(settings.DataFile,settings.DataThinning,settings.TargetChromosome,0.99);
	auto path = GetHarmonics(d,settings,gp);
	
	Data d2 = Data(settings.DataFile,settings.DataThinning,settings.TargetChromosome,settings.MemorySmoothing);
	auto path2 = GetHarmonics(d2,settings,gp);
	basicPlot(gp,d,0);
	TransitionPlot(gp,d,path,"Filtered");
	TransitionPlot(gp,d2,path2,"Raw");
	gp.SetLegend(true);
	gp.Show();

	Log("Deforest routine completed. Have a nice day.\n\n")
}