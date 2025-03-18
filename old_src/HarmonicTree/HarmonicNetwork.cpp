#include "HarmonicNetwork.h"

HarmonicNetwork::HarmonicNetwork(const Data & d, int assignedChromosome,const Settings & s,bool scan)
{
	Initialise(d,assignedChromosome,s,scan);
}

void HarmonicNetwork::Initialise(const Data & d, int assignedChromosome, const Settings & s,bool scan)
{
	Memory = s.L;
	Qmax = s.Qmax;
	ScanMode = scan;
	int stepSize = (d.Chromosomes[assignedChromosome].Idx[1] - d.Chromosomes[assignedChromosome].Idx[0]);
	BufferSize = s.L/stepSize+1; // +1 to allow for either-side grabbing
	
	
	
	DataSize = d.Chromosomes[assignedChromosome].Counts.size();
	MyChromosome = assignedChromosome;

	logPloidyPrior = log(s.PloidyPrior);
	logContinuityPrior = log(s.ContinuityPrior);
	Ploidy = s.Ploidy;
	Paths.resize(BufferSize);

	for (int i = 0; i < BufferSize; ++i)
	{
		Paths[i].resize(Qmax);
	}
}


// //some initial overload of the probability function
// double logP(int k, int q, double nu, double gamma)
// {
// 	double d = (k - q*nu)/gamma;

// 	return log(1.0/(M_PI * gamma * (1 + d*d)));
// 	// return -0.5 * d* d;
// }

Path HarmonicNetwork::Navigate(const Data & d, ProbabilityModel & prob)
{
	//do the initial layer first
	// 
	int k0 = d.Chromosomes[MyChromosome].Counts[0];
	

	for (int q = 0; q < Qmax; ++q)
	{
		double score = prob.logP(k0,q);
		if (q != Ploidy)
		{
			score += logPloidyPrior;
		}
		Paths[0][q].InitialStep(score,q);

	}
	
	for (chr_int i = 1; i < DataSize; ++i)
	{
		int bufferIndex = i % BufferSize;
		
		int lowerIndex = bufferIndex-1;
		if (lowerIndex < 0)
		{
			lowerIndex = BufferSize-1;
		}


		int k = d.Chromosomes[MyChromosome].Counts[i];
		//start by grabbing the simple continuity option
		for (int q = 0; q < Qmax; ++q)
		{
			double nodeProb = prob.logP(k,q);
			if (q!= Ploidy)
			{
				nodeProb += logPloidyPrior;
			}
			Path * BestPath = &Paths[lowerIndex][q];
			double BestScore = BestPath->Score + nodeProb;
			chr_int bestBack = -1;
			chr_int bestIndex = d.Chromosomes[MyChromosome].Idx[i];
			double cumProb = BestPath->CumulativeLinearScore + nodeProb;
			bool jumped=false;
			
			//do some additional checking
			if (i >= BufferSize-1)
			{
				int jumpIndex = bufferIndex + 1;
				if (jumpIndex >= BufferSize)
				{
					jumpIndex = 0;
				}
				
				for (int qq = 0; qq < Qmax; ++qq)
				{
					if (qq != q)
					{
						double cumulativeDifference = cumProb - Paths[jumpIndex][q].CumulativeLinearScore;


						double jumpScore = Paths[jumpIndex][qq].Score + cumulativeDifference + logContinuityPrior;

						if (jumpScore > BestScore)
						{
							BestScore = jumpScore;
							jumped =true;
							BestPath = &Paths[jumpIndex][qq];
							bestBack = Paths[jumpIndex][qq].Idx;
						}
					}
				}
			}
			if (ScanMode)
			{
				Paths[bufferIndex][q].RecordScore(BestScore,q,bestIndex);
			}
			else
			{
				Paths[bufferIndex][q].AddStep(BestScore,q,bestIndex,BestPath,bestBack);
			}
			Paths[bufferIndex][q].CumulativeLinearScore = cumProb; // ensure that the linear accumulation stays in the right place
		}	
	}

	return BestPath();
}

Path HarmonicNetwork::BestPath()
{
	// std::cout << "Attempting to find best path" << std::endl;
	int finalIdx = (DataSize-1) % BufferSize;

	double bestScore;
	Path * BestPath;
	double p = 0;
	for (int q = 0; q < Qmax; ++q)
	{
		double score = Paths[finalIdx][q].Score;
		p = score;
		if (q == 0 || score > bestScore)
		{
			bestScore = score;
			BestPath = &Paths[finalIdx][q];
		}
	}

	return *BestPath;
}