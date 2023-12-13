#include "GetHarmonics.h"



Path GetHarmonics(Data & d, Settings & s,JSL::gnuplot & gp)
{
	double gamma = 0.3;

	
	std::vector<HarmonicNetwork> Ns(d.Chromosomes.size());
	for (int c = 0; c < Ns.size(); ++c)
	{
		Ns[c].Initialise(d,c,s,true);
	}


	int count = 0;
	double bN = 30;
	double bestScore = -99999999999;
	// JSL::ProgressBar pb(100);
	
	Distributor dist(s.ParallelWorkers,d,Ns);
	for (double nus = d.Mean*0.4; nus < d.Mean*0.6; nus+=0.2*d.Mean/300)
	{
		dist.UpdateParameters(nus,gamma);
		dist.Signal(2);
		for (int i = 0; i < dist.MainChromAssigment.size(); ++i)
		{
			int chrom = dist.MainChromAssigment[i];
			Ns[chrom].Navigate(d,nus,gamma);
		}
		dist.Gather();
		// for (int c = 0; c < Ns.size(); ++c)
		// {
			
		// 	Path best = Ns[c].Navigate(d,nus,gamma);
		// 	if (c == 0)
		// 	{
		// 		if (best.Score > bestScore)
		// 		{
		// 			bN = nus;
		// 			bestScore = best.Score;
		// 		}
		// 	}
			
		// }


		// pb.Update(count);
	}
	Ns[0].ScanMode = false;
	Ns[0].Navigate(d,bN,gamma);
	Path best = Ns[0].BestPath();
	best.Nu = bN;	
	return best;
}

