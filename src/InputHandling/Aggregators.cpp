#include "Aggregators.h"

void SplitCheck(const std::string & line,int n)
{
	if (n != 3)
	{
		throw std::runtime_error("Error reading line '" + line + "'Line splits into" + std::to_string(line.size()) + " fields. Expected 3 fields corresponding to samtools depth output. You likely got the wrong delimiter");
	}
}

DataHolder AggregateStream(std::istream& inputStream)
{
	
	std::ostringstream chromosomeManifest("");
	std::vector<int> standardWindows = {100,1000,10000};
	if (JSL::FindXInY(Settings.AccumulationFactor,standardWindows) == -1)
	{
		standardWindows.push_back(Settings.AccumulationFactor);
	}
	std::vector<Aggregator> crawler;
	JAR::Archive tar;
	if (Settings.CreateArchive)
	{
		LOG(INFO) << "Creating archive during memory read process";
		std::string outname = Settings.Output + ".gca";
		tar.Open(outname,std::ios::out);
		for (auto window : standardWindows)
		{
			crawler.push_back(Aggregator(window,tar));
		}
	}
	else
	{
		LOG(INFO) << "Archive file not being created; reading stream into memory";
	}
	int accumulator = Settings.AccumulationFactor;
	char delim = Settings.StreamDelimiter;
	LOG(DEBUG) << "Stream delimiter is '" << delim << "', accumulation factor=" << accumulator;

	//initialise output vector
	DataHolder data;

	//various counters and trackers
	int count = 0; //the counterpart to accumulator
	int chr = -1; //index of current chromosome
	int cidx = 0; //data index within current chromosome (distinct from base index: cidx = idx % accumulator) 
	std::string previousChromosome = "";
	
	std::string PIPE_LINE;
	while (std::getline(inputStream,PIPE_LINE))
	{
		std::vector<std::string> line = JSL::split(PIPE_LINE,delim);

		SplitCheck(PIPE_LINE,line.size());

		dnaindex idx = std::stoll(line[1]);
		int k = std::stoi(line[2]);
		if (line[0] != previousChromosome)
		{
			previousChromosome = line[0];
			chromosomeManifest << previousChromosome << "\n";
			for (int i = 0; i < crawler.size(); ++i)
			{
				std::string name = previousChromosome + "_" + std::to_string(standardWindows[i]) + ".dat";
				crawler[i].NewFile(name);
			}
			data.push_back(CoverageArray());
			chr+=1;
			cidx = -1;
			count = 0;
			LOG(INFO) << "Scanning new chromosome " << line[0];
			
		}	
		if (count % accumulator == 0)
		{
			cidx += 1;
			data[chr].AddData(idx,k);
			count = 1;
		}
		else
		{
			data[chr][cidx].Coverage += k;
			data[chr][cidx].SquareSum += k*k;
			++count;
		}
		for (int i = 0; i < crawler.size(); ++i)
		{
			crawler[i].Update(idx,k);
		}
	}
	if (Settings.CreateArchive)
	{
		LOG(INFO) << "Writing manifest to file";
		tar.Write(MANIFEST_FILE_NAME,chromosomeManifest.str());
		for (int i = 0; i < crawler.size(); ++i)
		{
			crawler[i].Flush();
		}
	}
	return data;
}

