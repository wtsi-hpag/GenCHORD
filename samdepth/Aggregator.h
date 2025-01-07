#include <fstream>
#include "tar_to_stream.h"
struct Aggregator
{
	std::ofstream & Tar;
	std::string File;
	int Counter;
	int Sum;
	int HopSize;
	long long int Position;
	int BufferCount = 0;
	std::ostringstream buffer;

	Aggregator(int hop,std::ofstream & tarStream) : Tar(tarStream)
	{
		HopSize = hop;
	}
	void NewFile(std::string fileName)
	{
		if (BufferCount != 0)
		{
			Flush();
		}
		File = fileName;
		Counter = 0;
		Sum = 0;
		Position = 0;
	}
	
	void Close()
	{
		Flush();
	}
	void Update(int idx, int k)
	{
		if (Position == 0)
		{
			Position = idx;
		}
		Sum += k;
		Counter += 1;
		if (Counter == HopSize)
		{
			buffer << Position << " " << Sum << "\n";
			BufferCount +=1;
			if (BufferCount > 10000)
			{
				BufferCount = 0;
				Flush();
				buffer.clear();
			}
			Counter = 0;
			Sum = 0;
			Position = 0;
		}
	}
	void Flush()
	{
		std::string temp = buffer.str();
		tar_to_stream(
			Tar, {
				.filename{File},
				.data{std::as_bytes(std::span{temp})},
			}
			
		);
	}
};