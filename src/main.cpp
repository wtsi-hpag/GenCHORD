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
		auto hist = Data.Histogram();

		double s= 0;
		Model P(hist.size(),6,Settings.AccumulationFactor);
		for (int i = 0; i < 2000;++i)
		{
			P.Compute();
			s+= P.Score(hist);
		// LOG(DEBUG) << P.Score(hist);
		}
		LOG(DEBUG) << s;
	}
	catch (const std::exception& e) 
	{
		LOG(ERROR) << "CRITICAL\n\t" << e.what();
		return 1;
	}
	LOG(INFO) << "GenCHORD Complete, Goodbye";
	return 0;
}