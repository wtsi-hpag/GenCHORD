#include <iostream>

#include "JSL.h"
#include "Utility/Log.h"
#include "settings.h"
#include "InputHandling/ParseHandler.h"
#include "Probability/Model.h"

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

int main(int argc, char ** argv)
{

	std::ios_base::sync_with_stdio(false);
	std::cin.tie(nullptr);
	try
	{
		Settings.Configure(argc,argv);
		ConfigureLogging();
		WelcomeWagon();
 
		DataHolder Data = ParseData();
		LOG(DEBUG) << "Data received in main";
		Data.Analyse();

		Settings.AccumulationFactor = 1;
		Model P(Settings.AccumulationFactor*100,6,Settings.AccumulationFactor);
		P.Parameters.Nu = 10;
		P.Parameters.Variance = 120;
		P.Parameters.Epsilon = 1e-10;
		P.Parameters.Weight = {0,1,0,0,0,0};
		P.Parameters.Contamination = {0.1,0,0,0,0,0};
		// P.Kmax = Settings.AccumulationFactor * 100;
		// P.NHarmonic = 6;
		P.Compute();
		std::vector<double> ks = JSL::Vector::intspace(0,P.Kmax,1);
		std::vector<double> ps(ks.size());
		double q =0;
		for (int k = 0; k < ks.size(); ++k)
		{
			ps[k] = exp(P[k]);
			q += ps[k];
			// LOG(DEBUG) << k << " " << ps[k] << " " << q;
		}
		JSL::gnuplot gp;
		gp.Plot(ks,ps);
		gp.SetGrid(true);
		// gp.SetYLog(true);
		gp.Show();

	}
	catch (const std::exception& e) 
	{
		LOG(ERROR) << "CRITICAL\n\t" << e.what();
		return 1;
	}
	LOG(INFO) << "GenCHORD Complete, Goodbye";
	return 0;
}