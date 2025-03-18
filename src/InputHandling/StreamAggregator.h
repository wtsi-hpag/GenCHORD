#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include "Data.h"
#include "../settings.h"
#include "../Utility/Log.h"
#include "Archiver.h"

const char ARCHIVE_DELIMITER = ' ';


struct StreamAggregator
{
	JAR::Archive & Tar;
	std::string File;
	int Counter;
	lint Sum;
	lint SqSum;
	int HopSize;
	lint Size;
	long long int Position;
	int BufferCount = 0;
	std::ostringstream buffer;

	StreamAggregator(int hop,JAR::Archive & tarStream) : Tar(tarStream)
	{
		HopSize = hop;
	}
	void NewFile(std::string fileName)
	{
		if (BufferCount != 0)
		{
			Flush();
		}
		buffer.str("");
		File = fileName;
		Counter = 0;
		Sum = 0;
		SqSum = 0;
		Size = 0;
		Position = -1;	
	}
	
	void Update(int idx, int k)
	{
		
		if (Position == -1)
		{
			LOG(DEBUG) << "\t Crawling position " << idx;
			Position = idx;
		}
		Sum += k;
		SqSum += pow(k,2);
		if (SqSum < 0)
		{
			throw overflow_error("Aggregation Squaresum has overflown");
		}
		Counter += 1;
		if (Counter == HopSize)
		{
			++Size;
			buffer << Position << ARCHIVE_DELIMITER << Sum << ARCHIVE_DELIMITER << SqSum << "\n";
			
			double err = abs( pow((double)Sum,2) - (double)SqSum);

			BufferCount +=1;
			Counter = 0;
			Sum = 0;
			SqSum = 0;
			Position = -1;
		}
	}
	void Flush()
	{
		LOG(DEBUG) << "Flushing string of size " << buffer.str().size() << " (" << Size << " lines) to " << File;
		BufferCount = 0;
		Tar.Write(File,buffer.str());
	}
};