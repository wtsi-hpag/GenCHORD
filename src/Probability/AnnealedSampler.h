#pragma once
#include "Model.h"
#include "../settings.h"
#include "../InputHandling/Data.h"

class AnnealedSampler
{
	private:
		std::vector<int> Histogram;
		OptimiserParameters Vector;
		OptimiserParameters Proposed;
	public:
		AnnealedSampler(const DataHolder & data);

		Model Fit();
};