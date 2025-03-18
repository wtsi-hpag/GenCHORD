#define GNUPLOT_NO_TIDY 1
#include <iostream>
#include "JSL.h"
#include "Utility/Log.h"
#include "settings.h"
#include "InputHandling/ParseHandler.h"
// #include "Probability/Model.h"
#include "Probability/AnnealedSampler.h"
#include "Utility/Random.h"
#include "Utility/Timer.h"
#include "Harmonics/HarmonicTree.h"
void ConfigureLogging()
{
	std::ios_base::sync_with_stdio(false);

	LOGCFG.headers = Settings.LogHeaders;
	LOGCFG.SetLevel(Settings.LogLevel);

	
}

void WelcomeWagon()
{
	LOGCFG.headers = false;
	int n = 20;
	std::string title = " GenCHORD ";
	std::string subtitle = "Harmonic Genomic Rearrangement Analysis\n";
	int off = n - (int)((-title.size() + subtitle.size() -1)/2);
	LOG(INFO) << "\n" <<std::string(n,'=') << title <<  std::string(n,'=') << "\n" << std::string(off,' ') << subtitle;

	LOGCFG.headers = Settings.LogHeaders;

	if (Settings.LogLevel > 1)
	{
		LOG(INFO) << "Logging level set to " << Settings.LogLevel;
		LOG(ERROR) << "\tErrors appear like this";
		LOG(WARN) << "\tWarnings appear like this";
		LOG(INFO) << "\tInformation appears like this";
		LOG(DEBUG) << "\tDebug output appears like this";	
	}
}

void ProcessFunction()
{
	LOG(INFO) << "PROCESS MODE\n\tIn this mode, data is read in and processed, but no further action is taken.";
	if (!Settings.CreateArchive)
	{
		LOG(WARN) << "Archiving mode is deactivated. No Archive will be created.";
	}
	if (Settings.TreeMode)
	{
		LOG(WARN) << "PROCESS MODE and TREE MODE active. PROCESS MODE takes priority over TREE MODE";
		Settings.TreeMode = false;
	}
	DataHolder Data = ParseData();

	LOG(INFO) << "Parsing complete. This completes PROCESS MODE.";
	return;
}

void TreeFunction()
{
	LOG(INFO) << "TREE MODE\n\tThis is the standard mode for genchord. The data will be read in, a probability model fitted, and a GenTree generated";
	

	DataHolder Data = ParseData();
	Data.Analyse();

	AnnealedSampler AS(Data);
	
	auto model = AS.Fit();
	model.PrepareHarmonics();
	LOG(INFO) << "Harmonics prepared";
	int C = 6;
	JSL::gnuplot gp;
	int row = 2;
	int col = 3;
	gp.SetMultiplot(row,col);
	gp.WindowSize(2000,900);
	gp.SetSuperTitle("Accumulation Length: " + std::to_string(Settings.AccumulationFactor));
	// y = 0;
	// x = 0;
	// for (int c = 0; c < Data.size(); ++c)
	for (int c =0; c < C; ++c)
	{
		gp.SetAxis(c);
		gp.SetTitle("Chromosome " + Data[c].Name);
		LOG(INFO) << "Beginning fine-tuning of Chromosome " << Data[c].Name;
		auto tunedModel = AS.FineTune(model,c);

		LOG(INFO) << "I find the contamination to be " << tunedModel.Parameters.Eta;
		std::string s = "(";
		for (int i = 0; i < tunedModel.Parameters.Contamination.size(); ++i)
		{
			s += std::to_string(tunedModel.Parameters.Contamination[i]/(Settings.Ploidy - i)) + ", ";
		}
		s+=")";
		LOG(DEBUG) << "\tPer-node Contamination: " << s;
		HarmonicTree Network(tunedModel,Data,c);
		auto path = Network.Navigate();

		auto xy = Data[c].GetCoverage();
		std::vector<double> yRound;
		std::vector<double> xRound;
		int skipper = max(1,(int)xy.X.size()/10000);
		for (int i = 0; i < xy.Y.size(); ++i)
		{
			if (i % skipper == 0)
			{
				xRound.push_back(xy.X[i]);
				yRound.push_back(xy.Y[i] * 1.0/(tunedModel.Parameters.Nu * Settings.AccumulationFactor));
			}
		}
		// JSL::gnuplot gp;
		gp.Plot(xRound,yRound);
		gp.SetTitle("Chromosome " + Data[c].Name);
		gp.SetYRange(0,model.Kmax+10);

		std::vector<int> x= {0};
		std::vector<double> y = {0};
		double prevy = 0;
		double maxHarmonic = 0;
		for (auto pos : path.Route)
		{
			x.push_back(pos.Index);
			x.push_back(pos.Index);
			y.push_back(prevy);
			int q = pos.Value;
			prevy = (q + tunedModel.Parameters.Contamination[q]);
			if (prevy > maxHarmonic)
			{
				maxHarmonic = prevy;
			} 
			y.push_back(prevy);
		}
		x.push_back(Data[c][Data[c].size()-1].Index);
		y.push_back(prevy);
		gp.Plot(x,y);
		gp.SetYRange(0,int(maxHarmonic+1.5));
		gp.SetGrid(true);
		gp.SetXLabel("Index");
		gp.SetYLabel("Copy Number (Estimated)");
	}
	gp.Show();


}

int main(int argc, char ** argv)
{
	srand(0);
	std::ios_base::sync_with_stdio(false);
	std::cin.tie(nullptr);
	try
	{
		Settings.Configure(argc,argv);
		ConfigureLogging();
		WelcomeWagon();
 

		if (Settings.ProcessMode)
		{
			ProcessFunction();
		}

		if (Settings.TreeMode)
		{
			TreeFunction();
		}		
	}
	catch (const std::exception& e) 
	{
		LOG(ERROR) << "CRITICAL\n\t" << e.what();
		return 1;
	}
	LOG(INFO) << "GenCHORD Complete, Goodbye.";
	return 0;
}