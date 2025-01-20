#include "StreamAggregator.h"
#include "../settings.h"
#include "../Utility/Log.h"
#include "Archiver.h"

const std::string MANIFEST_FILE_NAME = "chromosome.manifest";


//the basic function which scans over the input stream and saves it to memory. No archiving or other fancy business
DataHolder ParseRawInput(std::istream & inputStream);


struct Counters //things which tick up to track current position & trigger behaviour at certain points
{
	int ChromosomeCount;
	int DataIndex;
	int SumCount;
	int IndependenceCount;
	Counters(){
		ChromosomeCount = -1;
		DataIndex = -1;
		SumCount = 0;
		IndependenceCount=0;
	}
	void ChromosomeReset()
	{
		++ChromosomeCount;
		DataIndex = -1;
		SumCount = 0;
		IndependenceCount = 0;
	}
};

class StatefulReader
{
	std::ostringstream chromosomeManifest;
	std::vector<StreamAggregator> & Crawler;
	std::vector<CoverageArray> & Data;
	JAR::Archive & Tar;
	Counters Count;
	std::string PreviousChromosome;
	bool ReadingValidChromosome;
	dnaindex PreviousIndex;
	int ExpectedGap;
	StatefulReader(std::vector<CoverageArray> & data, std::vector<StreamAggregator> & crawler, JAR::Archive tar,int gap) : Data(data), Crawler(crawler), Tar(tar)
	{
		Count = Counters();
		chromosomeManifest.str("");
		PreviousChromosome = "";
		ReadingValidChromosome = true;
		ExpectedGap = gap;
		PreviousIndex = -gap;
	}

	void CheckNewChromosome(std::string chromosome);
	void ChromosomeCleanup();

}