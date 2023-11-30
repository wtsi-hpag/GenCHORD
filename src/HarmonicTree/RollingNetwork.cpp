#include "RollingNetwork.h"

RollingNetwork::RollingNetwork(const Data & d, int assignedChromosome, int qMax, int L)
{
	Initialise(d,assignedChromosome,qMax,L);
}

void RollingNetwork::Initialise(const Data & d, int assignedChromosome, int qMax, int L)
{
	Memory = L;
	Qmax = qMax;
	int stepSize = (d.Chromosomes[assignedChromosome].Idx[1] - d.Chromosomes[assignedChromosome].Idx[0]);
	BufferSize = L/stepSize+1; // +1 to allow for either-side grabbing
	DataSize = d.Chromosomes[assignedChromosome].Counts.size();
	MyChromosome = assignedChromosome;
	Paths.resize(BufferSize);

	for (int i = 0; i < BufferSize; ++i)
	{
		Paths[i].resize(Qmax);
	}
}


//some initial overload of the probability function
double logP(int k, int q, double nu, double gamma)
{
	double d = (k - q*nu)/gamma;

	return -0.5 * d* d;
}

void RollingNetwork::Navigate(const Data & d, double nu, double gamma)
{

	//do the initial layer first
	// 
	int k0 = d.Chromosomes[MyChromosome].Counts[0];

	double logContinuityPrior = 0;
	double logPloidyPrior = log(1);

	//reset state to zero
	for (int i = 1; i < BufferSize; ++i)
	{
		for (int q = 0; q < Qmax; ++q)
		{
			Paths[i][q].InitialStep(-1,-1);
		}
	}
	for (int q = 0; q < Qmax; ++q)
	{
		double score = logP(k0,q,nu,gamma);
		if (q != 2) // set the ploidy settings better
		{
			score += logPloidyPrior;
		}
		// Paths[0][q].AddStep(score,q,0);
		// Paths[0][q].CumulativeLinearScore = score;
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
			double nodeProb = logP(k,q,nu,gamma);
			if (q!= 2)
			{
				nodeProb += logPloidyPrior;
			}
			Path * BestPath = &Paths[lowerIndex][q];
			double BestScore = BestPath->Score + nodeProb;
			chr_int bestBack = -1;
			chr_int bestIndex = d.Chromosomes[MyChromosome].Idx[i];
			// std::cout << "\t" << BestPath->Score << "  " << nodeProb << "  " << BestScore << std::endl;
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
							// std::cout << "Gonna jump " << std::endl;
						}

					}
				}

			}

			Paths[bufferIndex][q].AddStep(BestScore,q,bestIndex,BestPath,bestBack);
			Paths[bufferIndex][q].CumulativeLinearScore = cumProb; // ensure that the linear accumulation stays in the right place
		}

		
		// std::cout << i << "   " << std::setprecision(10) <<  Paths[bufferIndex][1].Score << "  " << Paths[bufferIndex][1].CumulativeLinearScore << std::endl;
	}


}

Path RollingNetwork::BestPath()
{
	// std::cout << "Attempting to find best path" << std::endl;
	int finalIdx = (DataSize-1) % BufferSize;

	double bestScore;
	Path * BestPath;
	double p = 0;
	for (int q = 0; q < Qmax; ++q)
	{
		double score = Paths[finalIdx][q].Score;

		// std::cout << "Path " << q << " has " << score;

		// if (q > 0)
		// {
		// 	std::cout << " (delta = " << score - p << ")";
		// }
		// std::cout << std::endl;
		p = score;
		if (q == 0 || score > bestScore)
		{
			bestScore = score;
			BestPath = &Paths[finalIdx][q];
		}
	}

	return *BestPath;
}