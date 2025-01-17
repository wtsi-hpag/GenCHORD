#include <iostream>

#include "JSL.h"
#include "Utility/Log.h"
#include "settings.h"
#include "InputHandling/Streamer.h"

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


		// for (int i = 0; i < Data.size(); ++i)
		// {
		// 	int skipFactor = max(1, Data[i].size()/5000);
		// 	LOG(WARN) << skipFactor;
		// 	auto cov = Data[i].GetCoverage(skipFactor);
		// 	// std::vector<int> x = JSL::Vector::intspace(0,cov.size()*skipFactor-1,skipFactor);
			
		// 	JSL::gnuplot gp;
		// 	gp.SetTitle(Data[i].Name);
		// 	gp.Plot(cov.X,cov.Y);
		// 	gp.SetYRange(0,Data[i].AggregateMean*2);
		// 	gp.Show();
		// }
		
		// auto vec = Data.Histogram();
		
		// std::vector<int> x = JSL::Vector::intspace(0,vec.size()-1,1);
		// LOG(DEBUG) << vec.size() << " " << x.size();
		// JSL::gnuplot gp;
		// gp.Plot(x,vec);
		// gp.SetYLog(true);
		// gp.SetXLog(true);
		// gp.SetYRange(0.5,1e3);
		// gp.Show();
		
	}
	catch (const std::exception& e) 
	{
		LOG(ERROR) << "CRITICAL\n\t" << e.what();
		return 1;
	}
	LOG(INFO) << "GenCHORD Complete, Goodbye";
	return 0;
}