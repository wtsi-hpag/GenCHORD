#include "GetHarmonics.h"

Path GetHarmonics(Data & d, Settings & s,JSL::gnuplot & gp)
{
	double gamma = 0.3;

	auto start = std::chrono::high_resolution_clock::now();
	std::vector<HarmonicNetwork> Ns(d.Chromosomes.size());
	for (int c = 0; c < Ns.size(); ++c)
	{
		Ns[c].Initialise(d,c,s,true);
	}


	int count = 0;
	double bN = -1;
	double bestScore = -99999999999;
	JSL::ProgressBar pb(100);
	for (double nus = d.Mean*0.4; nus < d.Mean*0.6; nus+=0.2*d.Mean/100)
	{
		++count;
		
		for (int c = 0; c < Ns.size(); ++c)
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
			}
			
		}
		pb.Update(count);
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Rolling = " << ((double)std::chrono::duration_cast<std::chrono::microseconds>(end-start).count())/(pow(10,6)*count) << std::endl;
	Ns[0].ScanMode = false;
	Ns[0].Navigate(d,bN,gamma);
	Path best = Ns[0].BestPath();
	best.Nu = bN;	
	return best;
}

