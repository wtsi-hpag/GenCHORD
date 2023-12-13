// #define GNUPLOT_NO_TIDY
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
	d = Data(settings.DataFile,settings.DataThinning,settings.TargetChromosome,settings.MemorySmoothing);
	gp.Plot(d.Chromosomes[0].Idx,d.Chromosomes[0].Counts);		
	auto path = GetHarmonics(d,settings,gp);

	TransitionPlot(gp,d,path,"Harmonic nu=" + std::to_string(path.Nu));
	
	gp.SetLegend(true);
	gp.Show();

	Log("Deforest routine completed. Have a nice day.\n\n")
}