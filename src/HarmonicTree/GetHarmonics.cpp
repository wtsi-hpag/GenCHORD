#include "GetHarmonics.h"



Path GetHarmonics(Data & d, Settings & s,JSL::gnuplot & gp)
{
	double gamma = 0.3;

	
	std::vector<HarmonicNetwork> Ns(d.Chromosomes.size());
	for (int c = 0; c < Ns.size(); ++c)
	{
		Ns[c].Initialise(d,c,s,true);
	}


	double bN = 30;
	double bestScore = -99999999999;

	const int N  = s.nuResolution;
	JSL::ProgressBar pb(N);
	
	Distributor dist(s.ParallelWorkers,d,Ns);
	int count = 0;
	for (double nus = d.Mean*0.4; nus < d.Mean*0.6; nus+=0.2*d.Mean/N)
	{
		++count;
		dist.UpdateParameters(nus,gamma);
		dist.Signal(2);
		for (int i = 0; i < dist.MainChromAssigment.size(); ++i)
		{
			int chrom = dist.MainChromAssigment[i];
			Ns[chrom].Navigate(d,nus,gamma);
		}
		dist.Gather();
		Path best = Ns[0].Navigate(d,nus,gamma);
		if (best.Score > bestScore)
		{
			bN = nus;
			bestScore = best.Score;
		}
		pb.Update(count);
	}
	Ns[0].ScanMode = false;
	Ns[0].Navigate(d,bN,gamma);
	Path best = Ns[0].BestPath();
	best.Nu = bN;	
	return best;
}

