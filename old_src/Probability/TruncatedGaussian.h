#pragma once
#include <vector>
#include "../Utility/basicFunctions.h"
class TruncatedGaussian
{
	private:
	public:
		int Truncator;
		std::vector<std::vector<double>> PureData;
		int KMax;
		TruncatedGaussian(int kMax, int bounder);
		void Populate(double gamma);
};