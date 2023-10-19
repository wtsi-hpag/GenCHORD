#pragma once
#include <vector>
#include <string>
#include "basicFunctions.h"
#include "JSL.h"

#define chr_int int //change this definition to long long int for very long chromsomed individuals

struct Chromosome
{
	std::string Name;
	std::vector<chr_int> Idx;
	std::vector<int> Counts;
	int maxIdx = 0;
	int maxK;
	Chromosome(){};
	Chromosome(std::string name){Name = name;}

	int GetMaxK(int maxVal)
	{
		int truncated = 0;
		for (int i = 0; i < Counts.size(); ++i)
		{
			if (Counts[i] > maxVal)
			{
				Counts[i] = maxVal;
				++truncated;
			}
			if (i == 0 || Counts[i] > maxK)
			{
				maxK = Counts[i];
			}
		}
		return truncated;
	}
	void Add(chr_int id, int count)
	{
		Idx.push_back(id);
		Counts.push_back(count);
		if (id > maxIdx)
		{
			maxIdx = id;
		}
	}
};


// void GapCheck(std::string dataFile)
// {
// 	long long int i = 0;
// 	std::vector<std::string> name;
// 	std::vector<chr_int> idx;
// 	std::vector<int> coverage;

// 	forLineVectorIn(target, ' ',
// 		std::string name = FILE_LINE_VECTOR[0];
// 		std::string 
// 	);
// }


struct Data
{
	double mean;
	std::vector<Chromosome> Chromosomes;
	int maxK=0;
	double Mean;
	double Deviation;
	Data(std::string target, int thinning, std::string targetChromosome)
	{
		Log("Loading data from " << target  << std::endl)
		
		long long int i = 0; //reasonable to assume total genome length exceeds INT_MAX, even if chr_int_max does not
		long long int j = 0;
		double accumulator = 0;
		double accumulatorSq = 0;
		int chromID = -1;
		bool detectingNormalGap = true;
		chr_int normalGap;
		chr_int prev = -1;
		int spoofCount = 0;
		long long int s = 0;
		std::string currentFlag = "";
		forLineVectorIn(target, ' ',
			++s;
			std::string chromFlag = FILE_LINE_VECTOR[0];
			if (chromFlag == targetChromosome || targetChromosome == "all" )
			{
				if (chromID == -1 || chromFlag != currentFlag)
				{
					i = 0; //reset the thinning counter on a new chrom - prevents weird offsets
					Chromosomes.push_back(Chromosome(chromFlag));
					currentFlag = chromFlag;
					Log("\tLoading " << chromFlag << " at file line " << s << "\n");
					++chromID;
				}

				//spoof in expected zeros if there are gaps in the coverage file
				chr_int id=std::stoi(FILE_LINE_VECTOR[1]);
				int gap = 0;
				if (detectingNormalGap == true)
				{
					if (prev > -1)
					{	
						normalGap = id - prev;
						detectingNormalGap = false;
					}
				}
				else
				{
					gap = id - prev;
					if (gap != normalGap)
					{
						int spoofer = prev+normalGap;
						while (spoofer < id)
						{
							if (i % thinning == 0 )
							{
								Chromosomes[chromID].Add(spoofer,0);
								++spoofCount;
								++j;
							}
							++i;
							spoofer += normalGap;
						}
					}
				}
				prev = id;

				if (i % thinning == 0 )
				{
					
					int count = std::stoi(FILE_LINE_VECTOR[2]);
					accumulator += count;
					accumulatorSq += count * count;
					Chromosomes[chromID].Add(id,count);
					++j;
					if (count > maxK)
					{
						maxK = count;
					}
				}
				++i;
			}
		
		);
		
		Mean = accumulator/j;
		Deviation = sqrt(accumulatorSq/j - Mean * Mean); 

		int LudicrousValue = Mean + 20*Deviation;
		int truncated = 0;
		for (int i = 0; i < Chromosomes.size(); ++i)
		{
			truncated += Chromosomes[i].GetMaxK(LudicrousValue);
		}
		Log("\tI am asserting that any coverage above " << LudicrousValue << " is spurious. " << truncated << " Datapoints were affected\n")
		Log("\tSpoofed " << spoofCount << " missing entries as gaps\n\tLoaded of " << j << " datapoints\n\tData has global coverage  "<< Mean << "±" << Deviation  << "\n\tMaximum coverage value is " << maxK << std::endl;);
		
	}
};