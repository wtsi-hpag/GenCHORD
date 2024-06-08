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

template<class T>
void WriteMeta(std::ostringstream &stream, std::string name,T val)
{
	stream << name << ": " << val << "\n";
}


int main(int argc, char**argv)
{
	//read the command line arguments first
	JSL::Argument<std::string> configFile("__none__","config",argc,argv);
	Settings settings(configFile,argc,argv);
	
	//some pretty introductory text
	globalVerbose = !settings.Quiet;
	Log("==========================================\n\tCoverage Deforesting\n==========================================" << std::endl);
	settings.DataThinning = 1e3;
	JSL::mkdir(settings.OutputDirectory);
	//load the data file -- either from file, or from a pipe
	Data d(settings);
	
	settings.DataThinning=1e5;
	Data d2(settings);
	Log("Preparing probability models:\n")
	int Kmax = 4*d.Mean + d.Deviation;
	int qmax = settings.Qmax;
	auto model = Models::NegativeBinomial();

	
	NormaliseModel(model,d,settings, Kmax,qmax);

	// model.SetSignalParameters(11.60,5.33);
	// model.SetNoiseParameters(10.1452,4.7,0.30);
	
	Kmax = d.maxK;
	model.SetDimensionality(Kmax,qmax);
	model.SetGrids();


	// std::vector<double> alphas = {1e-10};
	// std::vector<int> L = {(int)3e5};
	// std::vector<double> ploidy = {0.1};
	std::vector<double> alphas = {1e-15,1e-10,1e-5};
	std::vector<int> L = {(int)1e5,(int)3e5,(int)6e5,(int)2e6};
	std::vector<double> ploidy = {0.1,0.5};
	std::string orig = settings.OutputDirectory;
	for (int i = 0; i < alphas.size(); ++i)
	{
		for (int j = 0; j < L.size(); ++j)
		{
			for (int p = 0; p < ploidy.size(); ++p)
			{
				settings.ContinuityPrior = alphas[i];
				settings.L = L[j];
				settings.PloidyPrior = ploidy[p];

				Log("Set parameters to " << alphas[i] << " " << L[j] << " " << ploidy[p] << std::endl;)
				
				std::string outStr = "alpha" + std::to_string((int)abs(log10(alphas[i]))) + "_L" + std::to_string((int)(10*log10(L[j]))) + "_prior" + std::to_string((int)(100*ploidy[p]));
				settings.OutputDirectory = orig + "/" + outStr + "/";
				JSL::mkdir(settings.OutputDirectory);

				Log("Beginning network navigation" << std::endl;)
				auto paths = GetHarmonics(d,settings,model);

				Log("Navigation complete, writing to output\n")
				
				std::string treeFile= "";

				std::ostringstream metaData;
				WriteMeta(metaData,"Data_Origin",settings.DataFile);
				WriteMeta(metaData,"Chromosome_Mode",settings.TargetChromosome);
				WriteMeta(metaData,"Harmonic_Frequency",model.SignalMean);
				WriteMeta(metaData,"Harmonic_Deviation",model.SignalSigma);
				WriteMeta(metaData,"Noise_Frequency",model.NoiseMean);
				WriteMeta(metaData,"Noise_Deviation",model.NoiseSigma);
				WriteMeta(metaData,"Noise_Weight",model.NoiseWeight);
				WriteMeta(metaData,"Minimum_Jump",settings.L);
				WriteMeta(metaData,"Continuity_Prior",settings.ContinuityPrior);
				WriteMeta(metaData,"Ploidy_Prior",settings.PloidyPrior);

				treeFile += metaData.str();


				for (int i = 0; i < paths.size(); ++i)
				{
					std::string s = paths[i].TreeOutput(d.Chromosomes[i].Name);
					treeFile += s;
					Log("\tPlotting " << i + 1 <<std::endl;)

					JSL::gnuplot gp;
					basicPlot(gp,d2,i);
					gp.WindowSize(2000,1400);
					TransitionPlot(gp,d2.Chromosomes[i],paths[i],"Inferred Curve");
					gp.SetTerminal("pngcairo");
					gp.SetLegend(true);
					gp.SetOutput(settings.OutputDirectory + "/" + d.Chromosomes[i].Name + ".png");
					// gp.SetPersistence(true);
					gp.Show();
				}
				std::string treeFileName = settings.OutputDirectory + "/coverage.tree";
				JSL::initialiseFile(treeFileName);
				JSL::writeStringToFile(treeFileName,treeFile);

				Log("\tOutput saved to " << settings.OutputDirectory << "\n");
			}
		}
	}
	Log("Deforest routine completed.Have a nice day.\n\n")
}