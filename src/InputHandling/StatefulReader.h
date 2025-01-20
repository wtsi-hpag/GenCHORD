#include "StreamAggregator.h"
#include "../settings.h"
#include "../Utility/Log.h"
#include "Archiver.h"

const std::string MANIFEST_FILE_NAME = "chromosome.manifest";

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
	std::vector<StreamAggregator> Crawler;
	JAR::Archive Tar;
	DataHolder & Data;
	Counters Count;
	std::string PreviousChromosome;
	bool ReadingValidChromosome;
	dnaindex PreviousIndex;
	int ExpectedGap;
	public:
		StatefulReader(DataHolder & data);

		void AddLine(const std::string & chromosome, const dnaindex & idx, const lint & k);
		void ChromosomeCleanup();
		void WriteManifest();
	private:
		void CheckNewChromosome(const std::string & chromosome);
		void InitialiseArchive();
		void IsValid(std::string chromosome);
		void ProcessLine(const dnaindex & idx, const lint & k);
};