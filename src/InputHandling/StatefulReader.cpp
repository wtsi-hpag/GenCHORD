#include "StatefulReader.h"
unsigned int zero = 0; //useful to have because casting sometimes causes issues

StatefulReader::StatefulReader(DataHolder & data) : Data(data)
{
	Count = Counters();
	chromosomeManifest.str("");
	PreviousChromosome = "";
	ReadingValidChromosome = false;
	ExpectedGap = Settings.DataGap;
	PreviousIndex = -ExpectedGap;

	InitialiseArchive();
}

void StatefulReader::InitialiseArchive()
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
		Tar.Open(outname,std::ios::out);
		
		for (auto window : standardWindows)
		{
			Crawler.push_back(StreamAggregator(window,Tar));
		}
	}
	else
	{
		LOG(INFO) << "Archive file not being created; reading stream into memory";
	}
}

void StatefulReader::WriteManifest()
{
	if (Settings.CreateArchive)
	{
		LOG(INFO) << "Writing manifest to file";
		Tar.Write(MANIFEST_FILE_NAME,chromosomeManifest.str());
	}
}

void StatefulReader::ChromosomeCleanup()
{
	if (ReadingValidChromosome)
	{
		if (Count.SumCount != 0)
		{
			Data[Count.ChromosomeCount].FlagTruncated(); //remove underfull elements
		}
		for (int i = 0; i < Crawler.size(); ++i)
		{
			Crawler[i].Flush();
		}
	}
}

void StatefulReader::IsValid(std::string chromosome)
{
	ReadingValidChromosome = (chromosome.find(Settings.IgnoreChromosomeFlag) == std::string::npos);
}



void StatefulReader::CheckNewChromosome(const std::string & chromosome)
{
	if (chromosome != PreviousChromosome)
	{
		//cleanup old chromosome (if valid)
		
		ChromosomeCleanup();

		IsValid(chromosome);

		if (!ReadingValidChromosome)
		{
			LOG(WARN) << "Ignoring chromosome " << chromosome << " (contains flag " << Settings.IgnoreChromosomeFlag << ")";
		}
		else
		{
			Data.NewChromosome();
			chromosomeManifest << chromosome << "\n";
			for (auto & crawl : Crawler)
			{
				crawl.NewFile(chromosome + "_" + std::to_string(crawl.HopSize) + ".dat");
			}
			Count.ChromosomeReset();
			PreviousIndex = -ExpectedGap;
			LOG(INFO) << "Scanning new chromosome: " << chromosome;
		}
		PreviousChromosome = chromosome;
	}
}


void StatefulReader::AddLine(const std::string & chromosome, const dnaindex & idx, const lint & k)
{
	CheckNewChromosome(chromosome);
	if (ReadingValidChromosome)
	{
		//without the -a command, samtools skips outputting large blocks of '0'; this detects those gaps and reinserts them.
		PreviousIndex += ExpectedGap;
		while (PreviousIndex < idx)
		{
			ProcessLine(PreviousIndex,zero);
			PreviousIndex += ExpectedGap;
		}
		ProcessLine(idx,k);

	}
}

void StatefulReader::ProcessLine(const dnaindex & idx, const lint & k)
{
	if (Count.IndependenceCount == 0)
	{
		auto & target = Data[Count.ChromosomeCount];
		if (Count.SumCount == 0)
		{
			++Count.DataIndex;
			target.AddData(idx,k);
		}
		else
		{
			target[Count.DataIndex].Coverage += k;
			target[Count.DataIndex].SquareSum += pow(k,2);
		}
		Count.SumCount = (Count.SumCount + 1) % Settings.AccumulationFactor;
		for (auto & crawl : Crawler)
		{
			crawl.Update(idx,k);
		}
	}
	Count.IndependenceCount = (Count.IndependenceCount + 1) % Settings.AutocorrelationLength;
}