#pragma once
#include <vector>
#include "Path.h"
#include "../settings.h"
#include "../DataFrame/data.h"
#include "../Probability/ProbabilityModel.h"

class HarmonicNetwork
{
	public:
		HarmonicNetwork(){};
		HarmonicNetwork(const Data & d, int assignedChromosome, const Settings & settings, bool scan);

		void Initialise(const Data & d, int assignedChromosome, const Settings & settings, bool scan);



		Path Navigate(const Data & d, ProbabilityModel & prob);

		Path BestPath();
		bool ScanMode;
	private:
		std::vector<std::vector<Path>> Paths;

		int Memory;
		int Qmax;
		int BufferSize;
		chr_int DataSize;
		int MyChromosome;
		double logContinuityPrior;
		double logPloidyPrior;
		int Ploidy;

};