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
		const DataHolder & Data;
		void Optimise(Model & model, int Nsteps);
	public:
		AnnealedSampler(const DataHolder & data);


		Model Fit();
		Model FineTune(Model original,int chromosome);
};