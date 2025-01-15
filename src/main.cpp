#include <iostream>
#include "JSL.h"
#include "Utility/Log.h"
#include "settings.h"
#include "InputHandling/Streamer.h"

void ConfigureLogging()
{
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
}

int main(int argc, char ** argv)
{

	try
	{
		Settings.Configure(argc,argv);
		ConfigureLogging();
		WelcomeWagon();
 
		DataHolder Data = ParseData();
		LOG(DEBUG) << "Data received in main";
		Data.Analyse();
		auto vec = Data.Histogram();
		
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