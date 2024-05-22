#pragma once
#include <vector>
#include <string>

#define chr_int int //change this definition to long long int for very long chromosomed individuals


struct ChromosomeCoverage
{
	std::string Name;
	std::vector<chr_int> Idx;
	std::vector<int> Counts;
	int maxIdx = 0;
	int maxK;
	ChromosomeCoverage(){};
	ChromosomeCoverage(std::string name){Name = name;}

	int RemoveSpuriousData(int maxVal);
	void Add(chr_int id, int count);
};
