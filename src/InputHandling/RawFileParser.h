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
		ChromosomeCount = 0;
		DataIndex = 0;
		SumCount = 0;
		IndependenceCount=0;
	}
};