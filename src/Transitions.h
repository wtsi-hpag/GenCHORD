#pragma once
#include "data.h"
#include <vector>

struct ChromTransitions
{
	std::vector<chr_int> Index;
	std::vector<int> Resonance;
	std::vector<double> Frequency;
};

struct Transitions
{
	std::vector<ChromTransitions> List;
	double Nu;
	double Gamma;
	double Sigma;
	double Score;
	Transitions(){};
	Transitions(int n)
	{
		List.resize(n);
	}
	void Add(int chrom,chr_int index, double nu, double sigma,int q)
	{
		List[chrom].Index.push_back(index);
		List[chrom].Resonance.push_back(q);
		List[chrom].Frequency.push_back(q*nu);
		Nu = nu;
		Sigma = sigma;
	}
	// void Plot(JSL::gnuplot & gp,Data & d)
	// {
	// 	std::vector<std::vector<chr_int>> idx(List.size());
	// 	std::vector<std::vector<double>> freq(List.size());
	// 	for (int c = 0; c < List.size(); ++c)
	// 	{
	// 		chr_int maxIdx = d.Chromosomes[c].maxIdx;
	// 		for (int i =0; i < List[c].Index.size(); ++i)
	// 		{
	// 			if (i > 0)
	// 			{
					
	// 				idx[c].push_back(List[c].Index[i]);
	// 				freq[c].push_back(List[c].Frequency[i-1]);
	// 				// res.push_back(List[c].Resonance[i-1]);
	// 			}
	// 			idx[c].push_back(List[c].Index[i]);
	// 			freq[c].push_back(List[c].Frequency[i]);
	// 			// res.push_back(List[c].Resonance[i]);
	// 		}
	// 		idx[c].push_back(maxIdx);
	// 		freq[c].push_back(freq[c][freq[c].size()-1]);
	// 		// res.push_back(res[res.size()-1]);
	// 	}
	// 	TransitionPlot(idx,freq,gp);
	// }
};
