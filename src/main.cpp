#include <iostream>

#include "JSL.h"
#include "Utility/Log.h"
#include "settings.h"
#include "InputHandling/ParseHandler.h"
// #include "Probability/Model.h"
#include "Probability/AnnealedSampler.h"
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
	
	AS.Fit();
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