#pragma once
#include <vector>
#include "../Utility/logFactorial.h"
#include "../Utility/basicFunctions.h"
#include <iostream>
#include <algorithm>


class NegativeBinomial
{
	private:
		double Sigma;
		int Kmax;
		LogFactorial LogFac;
		int Resolution;
	public:
		double dMu;
		double muMax;
		std::vector<std::vector<double>> PureData;
		
		NegativeBinomial(int kMax, int resolution,double maxMu);
		void Populate(double sigma);
		// ~NegativeBinomial();
};

