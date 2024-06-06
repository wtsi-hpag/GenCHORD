#include "GetHarmonics.h"



std::vector<Path> GetHarmonics(Data & d, Settings & s, ProbabilityModel & prob)
{

	Log("\tInitialising networks" << std::endl;);
	std::vector<HarmonicNetwork> Ns(d.Chromosomes.size());
	
	std::vector<Path> out;
	Log("\tBeginning Traversal\n")
	for (int c = 0; c < Ns.size(); ++c)
	{
		Log("\t\tChromosome " << d.Chromosomes[c].Name << std::endl;)
		Ns[c].Initialise(d,c,s,true);

		Ns[c].ScanMode = false;
		Ns[c].Navigate(d,prob);
		Path best = Ns[c].BestPath();
		best.Nu = prob.SignalMean;	
		out.push_back(best);
	}
	
	
	
	return out;
}

