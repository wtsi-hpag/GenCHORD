#pragma once
#include "../data.h"
#include "../settings.h"
#include "network.h"
#include <vector>
#include <chrono>
#include "RollingNetwork.h"

class BasicPath
{
	public:
		double Score;
		chr_int Length;
		std::vector<int> Values;
		std::vector<chr_int> Index;
		
		int vindex;
		int prev;

		BasicPath(){};
		BasicPath(chr_int L)
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

		std::vector<BasicPath> Paths;
		std::vector<BasicPath> PreviousPaths;
		int QMax;
		double Nu;
		double Gamma;
		Tree(int qMax, chr_int L, double nu, double gamma)
		{
			QMax = qMax;
			Paths = std::vector<BasicPath>(qMax,L);
			PreviousPaths = std::vector<BasicPath>(qMax,L);
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

		BasicPath BestPath()
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


double nu = 30.5;
double gamma = 5;

void RollingAssign(Data & d, Settings & s,JSL::gnuplot & gp)
{

	std::cout << "Attempting a network assign" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	std::vector<RollingNetwork> Ns(d.Chromosomes.size());

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

void TreeAssign(Data & d, Settings & s,JSL::gnuplot & gp)
{

	// std::vector<Network> Ns(d.Chromosomes.size());


	auto start = std::chrono::high_resolution_clock::now();

	std::vector<Network> Ns(d.Chromosomes.size());
	for (int c = 0; c < Ns.size(); ++c)
	{
		Ns[c].Init(d,c,s.Qmax);
	}

	int count = 0;
	for (double nus = 1; nus < 20; nus+=2)
	{
		++count;
		for (int c = 0; c < d.Chromosomes.size(); ++c)
		{
			Network & N = Ns[c];
			//assign start of network
			Node startNode;
			startNode.Assign(-1,-1);
			double logContinuityPrior = -10;
			double logPloidyPrior = log(s.PloidyPrior);
			for (int q = 0; q < s.Qmax; ++q)
			{
				int k = d.Chromosomes[c].Counts[0];
				double dist = (q * nu - k)/gamma;
				double p = -0.5*dist * dist;
				if (q != s.Ploidy)
				{
					p += logPloidyPrior;
				}
				N.Nodes[0][q].Score = p;
				N.Nodes[0][q].NodeProb = p;
				N.Nodes[0][q].CumulativeNodeProb = p;
				N.Nodes[0][q].Prev = &startNode;

			
			}
			int step = (d.Chromosomes[c].Idx[1] - d.Chromosomes[c].Idx[0]);
			int L_equiv = s.L/step;

			for (int i = 1; i < N.Nodes.size(); ++i)
			{
				//obvious connection
				int chromID = N.Nodes[i][0].Idx;//allows node resolution to be different to data
				int k =d.Chromosomes[c].Counts[chromID];
				for (int q = 0; q < s.Qmax; ++q)
				{
					double dist = (q * nu - k)/gamma;
					double nodeProb = -0.5*dist*dist;
					if (q != s.Ploidy)
					{
						// std::cout << "noq" << std::endl;
						nodeProb += logPloidyPrior;
					}
					N.Nodes[i][q].NodeProb = nodeProb;
					N.Nodes[i][q].CumulativeNodeProb = N.Nodes[i-1][q].CumulativeNodeProb + nodeProb;
					double bestProb = N.Nodes[i-1][q].Score + nodeProb;
					double stayProb = bestProb;
					Node * bestConnect = &N.Nodes[i-1][q]; //default is to assume constant is best
					bool jumped=false;
					if (i >= L_equiv)
					{
						for (int qq = 0; qq < s.Qmax; ++qq)
						{
							if (qq != q)
							{
								double jumpProb = logContinuityPrior + N.Nodes[i-L_equiv][qq].Score;
								
								double qDiff = N.Nodes[i][q].CumulativeNodeProb - N.Nodes[i-L_equiv][q].CumulativeNodeProb;

								jumpProb += qDiff;

								if (jumpProb > bestProb)
								{
									jumped=true;
									bestProb = jumpProb;
									bestConnect = &N.Nodes[i-L_equiv][qq];
								}
							}
						}
					}

					N.Nodes[i][q].Prev = bestConnect;
					N.Nodes[i][q].Score = bestProb;
					
				}
			}
			


			// if (nus == 1 && c == 0)
			// {

			// // gp.Plot(d.Chromosomes[0].Idx,d.Chromosomes[0].Counts);
			// // for (int i = 0; i < )
			// auto node = N.GetBestPath();
			// std::vector<int> id;
			// std::vector<int> harmonic;
			// std::vector<double> vals;
			// int prev;
			// while (node->Q != -1)
			// {
			// 	prev = node->Q;
			// 	if (id.size() == 0 || harmonic[harmonic.size()-1] != node->Q )
			// 	{
			// 		id.push_back(d.Chromosomes[0].Idx[node->Idx]);
					
			// 		harmonic.push_back(node->Q);
			// 		vals.push_back(node->Q*nu);
			// 		// std::cout << node->Idx << "  " << d.Chromosomes[0].Idx[node->Idx] << "  " << node->Q << std::endl;
			// 	}
			// 	node = node->Prev;
			// }

			// // id.push_back(0);
			// // harmonic.push_back(prev);
			// // vals.push_back(prev*nu);

			// std::reverse(id.begin(),id.end());
			// std::reverse(harmonic.begin(),harmonic.end());
			// std::reverse(vals.begin(),vals.end());
			// // gp.Scatter(id,vals);

			// std::vector<int> plotIdx;
			// std::vector<double> plotVals;
			// double prevID =0;
			// double prevVal = 0;
			// for (int i =0 ; i < vals.size(); ++i)
			// {
			// 	plotIdx.push_back(prevID);
			// 	plotVals.push_back(prevVal);
			// 	prevVal= vals[i];
			// 	plotIdx.push_back(prevID);
			// 	plotVals.push_back(prevVal);
			// 	prevID = id[i];
			// 	plotIdx.push_back(prevID);
			// 	plotVals.push_back(prevVal);
			// }
			// // gp.Scatter(id,vals,JSL::LineProperties::ScatterType(JSL::Star));
			// gp.Plot(plotIdx,plotVals);
			// gp.SetTitle("Brute Force");
			// gp.SetGrid(true);
			// }
			// std::cout << id.size() << std::endl;
		}
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << ((double)std::chrono::duration_cast<std::chrono::microseconds>(end-start).count())/(pow(10,6)*count) << std::endl;
	

	// double nu = 30;
	// double gamma = 1;
	// for (int i = 0; i < )
	// while (true)
	// {
	// 	std::cout << Ns.size() << "  " << sizeof(Node) << "  " << sizeof(Network) << "  " << alignof(double) << " " << alignof(int) << "  " << alignof(short)<< " " << alignof(Node *)<< std::endl;
	// }
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
