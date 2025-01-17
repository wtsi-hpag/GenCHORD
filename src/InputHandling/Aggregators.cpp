#include "Aggregators.h"

void SplitCheck(const std::string & line,int n)
{
	if (n != 3)
	{
		throw std::runtime_error("Error reading line '" + line + "'Line splits into" + std::to_string(line.size()) + " fields. Expected 3 fields corresponding to samtools depth output. You likely got the wrong delimiter");
	}
}

void parseLine(const std::string& line, char delim, std::string& chromosome, dnaindex& idx, unsigned int& k) {
    size_t firstDelim = line.find(delim);
    if (firstDelim == std::string::npos) throw std::runtime_error("Malformed line: " + line);

    size_t secondDelim = line.find(delim, firstDelim + 1);
    if (secondDelim == std::string::npos) throw std::runtime_error("Malformed line: " + line);

    chromosome = line.substr(0, firstDelim);

    auto idxParseResult = std::from_chars(line.data() + firstDelim + 1, line.data() + secondDelim, idx);
    if (idxParseResult.ec != std::errc()) throw std::runtime_error("Failed to parse idx");

    auto kParseResult = std::from_chars(line.data() + secondDelim + 1, line.data() + line.size(), k);
    if (kParseResult.ec != std::errc()) throw std::runtime_error("Failed to parse k");
}

void addLine(CoverageArray & data, std::vector<Aggregator> & crawler, int & cidx, int & count, dnaindex & idx, unsigned int &k)
{
	if (count == 0)
	{
		++cidx;
		data.AddData(idx, k);
	} 
	else
	{
		data[cidx].Coverage += k;
		data[cidx].SquareSum += pow(k,2);
	}
	count = (count + 1) % Settings.AccumulationFactor;


	for (int i = 0; i < crawler.size(); ++i)
	{
		crawler[i].Update(idx,k);
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
	std::vector<CoverageArray> data;

	//various counters and trackers
	int count = 0; //the counterpart to accumulator
	int chr = -1; //index of current chromosome
	int cidx = 0; //data index within current chromosome (distinct from base index: cidx = idx % accumulator) 
	std::string previousChromosome = "";
	int gap = Settings.DataGap;
	dnaindex prevIndex = -gap;
	std::string PIPE_LINE;
	unsigned int zero = 0;
	while (std::getline(inputStream,PIPE_LINE))
	{

		std::string chromosome;
        dnaindex idx;
        unsigned int k;

        parseLine(PIPE_LINE, Settings.StreamDelimiter, chromosome, idx, k);
		if (chromosome != previousChromosome)
		{
			previousChromosome = chromosome;
			chromosomeManifest << previousChromosome << "\n";
			for (int i = 0; i < crawler.size(); ++i)
			{
				std::string name = previousChromosome + "_" + std::to_string(standardWindows[i]) + ".dat";
				crawler[i].NewFile(name);
			}
			if (count != 0)
			{
				data[chr].FlagTruncated(); //remove underfull elements
			}
			data.push_back(CoverageArray(chromosome));
			chr+=1;
			cidx = -1;
			prevIndex = -gap;
			
			LOG(INFO) << "Scanning new chromosome: " << chromosome;
			count = 0;
			
		}	


		prevIndex+=gap;
		while (prevIndex < idx)
		{
			addLine(data[chr],crawler,cidx,count,prevIndex,zero);
			prevIndex+=gap;
		}
		addLine(data[chr],crawler,cidx,count,idx,k);

		
	}
	if (count != 0)
	{
		data[chr].FlagTruncated(); //remove underfull elements
	}
	if (Settings.CreateArchive)
	{
		
		for (int i = 0; i < crawler.size(); ++i)
		{
			crawler[i].Flush();
		}
		LOG(INFO) << "Writing manifest to file";
		tar.Write(MANIFEST_FILE_NAME,chromosomeManifest.str());
	}
	return data;
}

