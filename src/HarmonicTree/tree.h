#pragma once
#include "../data.h"
#include "../settings.h"
#include <vector>
#include <chrono>
#include "HarmonicNetwork.h"



double nu = 30.5;
double gamma = 5;

void RollingAssign(Data & d, Settings & s,JSL::gnuplot & gp)
{

	std::cout << "Attempting a network assign" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	std::vector<HarmonicNetwork> Ns(d.Chromosomes.size());

	for (int c = 0; c < Ns.size(); ++c)
	{
		Ns[c].Initialise(d,c,s.Qmax,s.L,true);
	}
	int count = 0;
	double bN = -1;
	double bestScore = -99999999999;
	for (double nus = 25; nus < 36; nus+=5)
	{
		++count;
	
		for (int c = 0; c < d.Chromosomes.size(); ++c)
		{
			Ns[c].Navigate(d,nus,gamma);
			Path best = Ns[c].BestPath();
			if (c == 0)
			{
				if (best.Score > bestScore)
				{
					bN = nus;
					bestScore = best.Score;
				}
			
				// std::cout << nus << "  " << best.Score << std::endl;
			}
		}
	}
	
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Rolling = " << ((double)std::chrono::duration_cast<std::chrono::microseconds>(end-start).count())/(pow(10,6)*count) << std::endl;
	std::cout << "\tBest nu = " << bN << std::endl;
	Ns[0].ScanMode = false;
	Ns[0].Navigate(d,bN,gamma);
	Path best = Ns[0].BestPath();

	
	std::vector<int> plotIdx = {best.Route[0].Index};
	double prev = best.Route[0].Value * nu;
	std::vector<double> plotVals = {prev};
	for (int i =1; i < best.Route.size(); ++i)
	{
		plotIdx.push_back(best.Route[i].Index);
		plotVals.push_back(prev);
		prev = best.Route[i].Value*nu;
		plotIdx.push_back(best.Route[i].Index);
		plotVals.push_back(prev);
	}

	// std::cout << JSL::Vector(best.PreviousIdx) << std::endl;
	// std::cout << JSL::Vector(best.Previous) << std::endl;
	plotIdx.push_back(d.Chromosomes[0].maxIdx);
	plotVals.push_back(prev);

	gp.Plot(d.Chromosomes[0].Idx,d.Chromosomes[0].Counts);

	gp.Plot(plotIdx,plotVals);
	// gp.Show();
}

