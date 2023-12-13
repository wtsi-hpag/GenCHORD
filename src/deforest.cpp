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

	std::vector<int> ws;
	std::vector<double> ts;
	for (int w = 0; w < 10; ++w)
	{
		settings.ParallelWorkers = w;
		auto start = std::chrono::high_resolution_clock::now();
		auto path = GetHarmonics(d,settings,gp);
		auto end = std::chrono::high_resolution_clock::now();
		double length = ((double)std::chrono::duration_cast<std::chrono::microseconds>(end-start).count())/(pow(10,6));
		TransitionPlot(gp,d,path,"Harmonic nu=" + std::to_string(path.Nu) + " w = " + std::to_string(w));
		ws.push_back(w);
		ts.push_back(length);
	}
	gp.SetLegend(true);
	gp.Show();

	JSL::gnuplot gp2;
	
	std::vector<double> theory = 1+JSL::Vector(ws);
	for (int i = 0; i < theory.size(); ++i)
	{
		theory[i] = ts[0] / theory[i];
	}
	gp2.Plot(ws,theory,JSL::LineProperties::Legend("Theoretical Max"));
	gp2.Plot(ws,ts,JSL::LineProperties::Legend("Real"));
	gp2.SetLegend(true);
	gp2.SetXLabel("Parallel Workers");
	gp2.SetYLabel("Execution Time");
	gp2.Show();
	
	Log("Deforest routine completed. Have a nice day.\n\n")
}