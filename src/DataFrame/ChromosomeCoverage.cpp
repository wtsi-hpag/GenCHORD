#include "ChromosomeCoverage.h"
#include <iostream>
int ChromosomeCoverage::RemoveSpuriousData(int maxVal)
{
	int truncated = 0;
	for (int i = 0; i < Counts.size(); ++i)
	{
		if (Counts[i] > maxVal)
		{
			Counts[i] = Counts[i-1];
			++truncated;
		}
		if (i == 0 || Counts[i] > maxK)
		{
			maxK = Counts[i];
		}
	}
	return truncated;
}

void ChromosomeCoverage::Add(chr_int id, int count)
{
	Idx.push_back(id);
	Counts.push_back(count);
	
	if (id > maxIdx)
	{
		maxIdx = id;
	}
}
void ChromosomeCoverage::Add(chr_int id, int count,double smooth)
{
	Idx.push_back(id);
	Counts.push_back(count);
	SmoothCounts.push_back(smooth);
	if (id > maxIdx)
	{
		maxIdx = id;
	}
}