#include "RawFileParser.h"
unsigned int zero = 0; //useful to have because casting sometimes causes issues

void parseLine(const std::string& line, char delim, std::string& chromosome, dnaindex& idx, unsigned int& k) 
{
	//check the file delimiters
    size_t firstDelim = line.find(delim);
    if (firstDelim == std::string::npos)
	{
		throw std::runtime_error("Malformed line: " + line);
	}

    size_t secondDelim = line.find(delim, firstDelim + 1);
    if (secondDelim == std::string::npos)
	{
		throw std::runtime_error("Malformed line: " + line);
	}

	//parse the results into the containers (throwing errors if broken)
    chromosome = line.substr(0, firstDelim);
    auto idxParseResult = std::from_chars(line.data() + firstDelim + 1, line.data() + secondDelim, idx);
    if (idxParseResult.ec != std::errc())
	{
		throw std::runtime_error("Failed to parse idx");
	}
    auto kParseResult = std::from_chars(line.data() + secondDelim + 1, line.data() + line.size(), k);
    if (kParseResult.ec != std::errc()) 
	{
		throw std::runtime_error("Failed to parse k");
	}
}

void addLine(CoverageArray & data, std::vector<StreamAggregator> & crawler, Counters & count, dnaindex & idx, unsigned int &k)
{
	if (count.SumCount== 0)
	{
		++count.DataIndex;
		data.AddData(idx, k);
	} 
	else
	{
		data[count.DataIndex].Coverage += k;
		data[count.DataIndex].SquareSum += pow(k,2);
	}
	count.SumCount = (count.SumCount + 1) % Settings.AccumulationFactor;

	for (int i = 0; i < crawler.size(); ++i)
	{
		crawler[i].Update(idx,k);
	}
}

void InitialiseArchive(std::vector<StreamAggregator> & crawler, JAR::Archive & tar)
{
	if (Settings.CreateArchive)
	{
		LOG(INFO) << "Creating archive during memory read process";
		
		std::vector<int> standardWindows = {100,1000,10000};
		if (JSL::FindXInY(Settings.AccumulationFactor,standardWindows) == -1)
		{
			standardWindows.push_back(Settings.AccumulationFactor);
		}
		std::string outname = Settings.Output + ".gca";
		tar.Open(outname,std::ios::out);
		
		for (auto window : standardWindows)
		{
			crawler.push_back(StreamAggregator(window,tar));
			
		}
		
	}
	else
	{
		LOG(INFO) << "Archive file not being created; reading stream into memory";
	}
}

void StatefulReader::ChromosomeCleanup()
{
	if (Count.SumCount != 0)
	{
		Data[count.ChromosomeCount].FlagTruncated(); //remove underfull elements
	}
	for (int i = 0; i < Crawler.size(); ++i)
	{
		Crawler[i].Flush();
	}
}

void StatefulReader::CheckNewChromosome(std::string chromosome)
{
	return (chromosome.find(Settings.IgnoreChromosomeFlag) == std::string::npos);
}




DataHolder ParseRawInput(std::istream& inputStream)
{
	std::ostringstream chromosomeManifest("");
	
	std::vector<StreamAggregator> crawler;
	JAR::Archive tar;
	InitialiseArchive(crawler,tar);
	std::vector<CoverageArray> data;

	StatefulReader Reader(data,crawler,tar,Settings.DataGap);


	int accumulator = Settings.AccumulationFactor;
	char delim = Settings.StreamDelimiter;
	LOG(DEBUG) << "Stream delimiter is '" << delim << "', accumulation factor=" << accumulator;

	//initialise output vector

	Counters count;

	std::string previousChromosome = "";
	int gap = Settings.DataGap;
	dnaindex prevIndex = -gap;
	
	bool validChromosome = false;
	
	std::string PIPE_LINE;
	while (std::getline(inputStream,PIPE_LINE))
	{
		std::string chromosome;
        dnaindex idx;
        unsigned int k;

        parseLine(PIPE_LINE, Settings.StreamDelimiter, chromosome, idx, k);

		//SNIP HERE
		if (chromosome != previousChromosome)
		{
			//clean up old chromosome
			if (validChromosome)
			{
				ChromosomeCleanup(data,crawler,count);	
			}

			
			//check & initialise a new chromosome entity
			validChromosome = IsValidChromosome(chromosome);
			if (!validChromosome)
			{
				LOG(WARN) << "Ignoring chromosome " << chromosome << " (contains flag " << Settings.IgnoreChromosomeFlag << ")";
			}
			else
			{
				data.push_back(CoverageArray(chromosome));

				chromosomeManifest << previousChromosome << "\n";
				for (int i = 0; i < crawler.size(); ++i)
				{
					std::string name = chromosome + "_" + std::to_string(crawler[i].HopSize) + ".dat";
					crawler[i].NewFile(name);
				}
			
				count.ChromosomeReset();
				prevIndex = -gap;
				LOG(INFO) << "Scanning new chromosome: " << chromosome;
			}

			previousChromosome = chromosome;
			
		}	
		if (validChromosome)
		{
			prevIndex+=gap;
			while (prevIndex < idx)
			{
				addLine(data[count.ChromosomeCount],crawler,count,prevIndex,zero);
				prevIndex+=gap;
			}
			addLine(data[count.ChromosomeCount],crawler,count,idx,k);
		}
	}
	//to here

	//finish off the currently open chromosome
	if (validChromosome)
	{
		if (count.SumCount != 0)
		{
			data[count.ChromosomeCount].FlagTruncated(); //remove underfull elements
		}

			
		for (int i = 0; i < crawler.size(); ++i)
		{
			crawler[i].Flush();
		}
	}

	//finalise the archive creation stuff
	if (Settings.CreateArchive)
	{
		LOG(INFO) << "Writing manifest to file";
		tar.Write(MANIFEST_FILE_NAME,chromosomeManifest.str());
	}
	return data;
}

