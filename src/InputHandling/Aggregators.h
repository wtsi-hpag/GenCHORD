#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include "Archiver.h"
#include "Data.h"
#include "../settings.h"
#include "../Utility/Log.h"
#include "Archiver.h"
//the basic function which scans over the input stream and saves it to memory. No archiving or other fancy business
std::vector<CoverageArray> AggregateStream(std::istream & inputStream);


struct Aggregator
{
	JAR::Archive & Tar;
	std::string File;
	int Counter;
	int Sum;
	int SqSum;
	int HopSize;
	long long int Position;
	int BufferCount = 0;
	std::ostringstream buffer;

	Aggregator(int hop,JAR::Archive & tarStream) : Tar(tarStream)
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
		Position = 0;	
	}
	
	void Update(int idx, int k)
	{
		if (Position == 0)
		{
			Position = idx;
		}
		Sum += k;
		SqSum += k*k;
		Counter += 1;
		if (Counter == HopSize)
		{
			buffer << Position << " " << Sum << " " << SqSum << "\n";
			BufferCount +=1;
			Counter = 0;
			Sum = 0;
			SqSum = 0;
			Position = 0;
		}
	}
	void Flush()
	{
		BufferCount = 0;
		Tar.Write(File,buffer.str());
	}
};