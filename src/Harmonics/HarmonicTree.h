#pragma once
#include "../settings.h"
#include "../Probability/Model.h"
#include "../InputHandling/Data.h"
#include "Path.h"
class HarmonicTree
{
	public:
		HarmonicTree(Model & model, DataHolder & Data, int chromosome);
		Path Navigate();

		Path BestPath();
	private:
		Model & Probability;
		DataHolder & Data;
		int TargetChromosome;

		

		std::vector<std::vector<Path>> Paths;

		int Memory;
		int Qmax;
		int BufferSize;
		dnaindex DataSize;
		double logContinuityPrior;
		double logPloidyPrior;
		int Ploidy;
};