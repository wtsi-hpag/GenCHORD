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
		std::vector<int> coverageArray(10,0);
		int totalData = 0;
		
		for (int i = 0; i < Data.size(); ++i)
		{
			int chromSize = Data[i].Size();
			totalData += chromSize;
			for (int j = 0; j < chromSize; ++j)
			{
				int k = Data[i][j].Coverage;
				
				if (k < 1000 * Settings.AccumulationFactor)
				{
					if (k >= coverageArray.size())
					{
						coverageArray.resize(k+10,0);
					}
					coverageArray[k] +=1;
				}
			}
		}
		double cumProb = 0.0;
		int cutOff = 0;
		while (cumProb < 0.999)
		{
			cumProb += coverageArray[cutOff] * 1.0/totalData;
			cutOff += 1;
			LOG(DEBUG) << cutOff << " " <<cumProb << " " << coverageArray.size();
		}
		coverageArray.resize(cutOff);
		JSL::gnuplot gp;
		double mod = 1.0;
		std::vector<int> v = JSL::Vector::linspace(0,(coverageArray.size()-1) * mod,coverageArray.size());
		gp.Plot(v,coverageArray);
		// gp.SetXLog(true);
		// gp.SetXRange(3000,5000);
		gp.SetYLog(true);
		gp.Show();
		// LOG(DEBUG) << coverageArray.size() << " is max coverage";
	}
	catch (const std::exception& e) 
	{
		LOG(ERROR) << "CRITICAL\n\t" << e.what();
		return 1;
	}
	LOG(INFO) << "GenCHORD Complete, Goodbye";
	return 0;
}