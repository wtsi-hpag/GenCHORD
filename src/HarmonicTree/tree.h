#pragma once
#include "../data.h"
#include "../settings.h"
#include "network.h"
#include <vector>
class Path
{
	public:
		double Score;
		chr_int Length;
		std::vector<int> Values;
		std::vector<chr_int> Index;
		
		int vindex;
		int prev;

		Path(){};
		Path(chr_int L)
		{
			Score = 0;
			Length = L;
			vindex = 0;
			prev = -1;
		}

		void Add(chr_int index,int value,double score)
		{
			Score += score;
			if (prev != value)
			{
				// std::cout << "New insertion " << prev << " -> " << value << " at " << index << " with distance " << TransitionDistance(index) << std::endl;
				Values.push_back(value);
				Index.push_back(index);
				++vindex;
			}
			prev = value;
		}

		int TransitionDistance(int index)
		{
			if (vindex < 2)
			{
				return INT_MAX;
			}
			else
			{
				return abs(index - Index[Index.size()-1]);
			}
		}
};

class Tree
{
	private:
		std::vector<int> q_buffer;
		std::vector<double> scoreBuffer;
	public:

		std::vector<Path> Paths;
		std::vector<Path> PreviousPaths;
		int QMax;
		double Nu;
		double Gamma;
		Tree(int qMax, chr_int L, double nu, double gamma)
		{
			QMax = qMax;
			Paths = std::vector<Path>(qMax,L);
			PreviousPaths = std::vector<Path>(qMax,L);
			Nu = nu;
			Gamma = gamma;
		}

		void AssignStep(int coord,int k)
		{
			PreviousPaths = Paths;
			
			for (int q = 0; q < QMax; ++q)
			{
				//work out best assignment
				int bestQ = -1;
				int L = 50000;
				double bestScore;

				

				bool forceFirst = true;
				
				for (int testQ = 0; testQ < QMax; ++testQ)
				{
					int distSinceTransition = Paths[testQ].TransitionDistance(coord);

					if (distSinceTransition > L || testQ == q)
					{
						
						double d = (k - testQ * Nu)/Gamma;
						double testScore = -0.5 * d*d;

						if (testQ != q)
						{
							testScore += -1;
						}
						if (testQ != 2)
						{
							testScore += log(0.4);
						}
						// if ()
						if (forceFirst || testScore > bestScore)
						{
							forceFirst = false;
							bestQ = testQ;
							bestScore = testScore;
							// std::cout << " New best test w " << distSinceTransition << "  " << testQ << "  " << q <<std::endl;
						}
					}
				}


				// std::cout << "Pushing " << bestQ << " to " << q << std::endl;
				if (bestQ!= q)
				{
					Paths[q] = PreviousPaths[bestQ];
				}
				Paths[q].Add(coord,q,bestScore);
			}
		}

		Path BestPath()
		{
			double bestScore = -99999;
			int bestIdx = 0;
			for (int i =0; i < Paths.size(); ++i)
			{
				if (i == 0 || Paths[i].Score > bestScore)
				{
					bestIdx = i;
					bestScore = Paths[i].Score;
				}
			}
			std::cout << "The path ending on " << bestIdx << " is best. \n\n Number of transitions: " << Paths[bestIdx].Values.size()  << std::endl;
			return Paths[bestIdx];
		}
};



void TreeAssign(Data & d, Settings & s)
{

	std::vector<Network> Ns(d.Chromosomes.size());
	for (int c = 0; c < d.Chromosomes.size(); ++c)
	{
		Ns[c] = Network(d,c,s.Qmax);
	}

	while (true)
	{
		std::cout << Ns.size() << "  " << sizeof(Node) << "  " << sizeof(Network) << "  " << alignof(double) << " " << alignof(int) << "  " << alignof(short)<< " " << alignof(Node *)<< std::endl;
	}
	// double nu = 30;
	// double gamma = 1;
	// int cm = 1;//d.Chromosomes.size();
	// for (int c = 0; c < cm; ++c)
	// {
	// 	std::cout << "Assigning chrom " << c+1 << std::endl;

	// 	int chromLength = d.Chromosomes[c].Counts.size();
	// 	Tree tree(s.Qmax,chromLength,nu,gamma);
	// 	for (int i = chromLength -1; i >= 0; --i)
	// 	{
	// 		tree.AssignStep(d.Chromosomes[c].Idx[i],d.Chromosomes[c].Counts[i]);
	// 		// std::cout << i << std::endl;
	// 	}

	// 	auto bp = tree.BestPath();



	// 	std::vector<int> plotIdx;
	// 	std::vector<double> plotVal;
	// 	double prev = 0;
	// 	int cs = 0;
	// 	for (int i = 0; i <bp.Index.size();++i)
	// 	{
	// 		plotIdx.push_back(bp.Index[i]);
	// 		plotIdx.push_back(bp.Index[i]);

	// 		plotVal.push_back(prev);
	// 		// std::cout << "Transition " << bp.Values[i] << "->" << (int)round(prev/nu) << " at " << bp.Index[i] << "(delta = " << abs(bp.Index[i] - cs) << ")"<< std::endl;
	// 		cs= bp.Index[i];
	// 		prev = bp.Values[i] * nu;
	// 		plotVal.push_back(prev);

	// 	}
	// 	plotIdx.push_back(0);
	// 	plotVal.push_back(prev);

	// 	JSL::gnuplot gp;
	// 	gp.Plot(d.Chromosomes[c].Idx,d.Chromosomes[c].Counts);
	// 	gp.Plot(plotIdx,plotVal);
	// 	gp.Show();
	// }
}
