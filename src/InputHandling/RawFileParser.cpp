#include "RawFileParser.h"


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

DataHolder ParseRawInput(std::istream& inputStream)
{
	DataHolder data;
	StatefulReader Reader(data);

	std::string PIPE_LINE;
	while (std::getline(inputStream,PIPE_LINE))
	{
		std::string chromosome;
        dnaindex idx;
        unsigned int k;

        parseLine(PIPE_LINE, Settings.StreamDelimiter, chromosome, idx, k);

		Reader.AddLine(chromosome,idx,k);
	}

	Reader.ChromosomeCleanup();
	Reader.WriteManifest();


	return data;
}

