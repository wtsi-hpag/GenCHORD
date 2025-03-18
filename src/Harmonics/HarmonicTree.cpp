#include "HarmonicTree.h"

HarmonicTree::HarmonicTree(Model & model, DataHolder & data, int chromosome) : Probability(model), Data(data)
{
	TargetChromosome = chromosome;

	Memory = Settings.MinJump;
	Qmax = Settings.HarmonicCount;
	int stepSize = (Data[chromosome][1].Index - Data[chromosome][0].Index);
	BufferSize = Memory/stepSize+1; // +1 to allow for either-side grabbing
	
	
	
	DataSize = Data[chromosome].size();

	logPloidyPrior = 0;//log(1.0 - model.Parameters.LogWeight[Settings.Ploidy]);
	logContinuityPrior = 0;//log(s.ContinuityPrior);
	Ploidy = Settings.Ploidy;
	Paths.resize(BufferSize);

	for (int i = 0; i < BufferSize; ++i)
	{
		Paths[i].resize(Qmax);
	}
}


Path HarmonicTree::Navigate()
{
	//do the initial layer first
	// 
	LOG(INFO) << "Beginning navigating chromosome " << Data[TargetChromosome].Name;
	int k0 = Data[TargetChromosome][0].Coverage;
	

	for (int q = 0; q < Qmax; ++q)
	{
		double score = Probability.HarmonicProbability(k0,q);
		if (q != Ploidy)
		{
			score += logPloidyPrior;
		}
		// score += Probability.Parameters.LogWeight[q];
		Paths[0][q].InitialStep(score,q);

	}
	LOG(DEBUG) << "\tInitial nodes set up";
	for (dnaindex i = 1; i < DataSize; ++i)
	{
		if (i < 10)
		{
			LOG(DEBUG) << i << " " << Data[TargetChromosome][i].Index;
		}
		int bufferIndex = i % BufferSize;
		
		int lowerIndex = bufferIndex-1;
		if (lowerIndex < 0)
		{
			lowerIndex = BufferSize-1;
		}


		int k = Data[TargetChromosome][i].Coverage;
		//start by grabbing the simple continuity option
		for (int q = 0; q < Qmax; ++q)
		{
			double nodeProb = Probability.HarmonicProbability(k,q);
			if (q!= Ploidy)
			{
				nodeProb += logPloidyPrior;
			}
			Path * BestPath = &Paths[lowerIndex][q];
			double BestScore = BestPath->Score + nodeProb;
			dnaindex bestBack = -1;
			dnaindex bestIndex = Data[TargetChromosome][i].Index;
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
			Paths[bufferIndex][q].AddStep(BestScore,q,bestIndex,BestPath,bestBack);
			
			Paths[bufferIndex][q].CumulativeLinearScore = cumProb; // ensure that the linear accumulation stays in the right place
		}	
	}

	return BestPath();
}

Path HarmonicTree::BestPath()
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