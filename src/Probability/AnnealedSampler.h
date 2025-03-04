#pragma once
#include "Model.h"
#include "../settings.h"
#include "../InputHandling/Data.h"

class AnnealedSampler
{
	private:
		std::vector<int> Histogram;
		OptimiserPack Vector;
		StateVector Proposed;

		// void Optimise(const DataHolder & data);
	public:
		AnnealedSampler(const DataHolder & data);


		Model Fit();
};