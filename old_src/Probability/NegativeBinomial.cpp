#include "NegativeBinomial.h"

NegativeBinomial::NegativeBinomial(int kMax, int resolution,double maxMu)
{
	Kmax = kMax;
	PureData.resize(kMax+1,std::vector<double>(resolution,0.0));
	dMu = maxMu * 1.0/resolution;
	Resolution = resolution;
}
void NegativeBinomial::Populate(double sigma)
{
	// sigma = 7;
	if (sigma < 1e-7)
	{
		sigma = 1e-7;
	}
	
	Sigma = sigma;

	
	double var = sigma*sigma;
	for (int k = 0; k <= Kmax; ++k)
	{
		//do the zeroth component because it involves log(0) if done naively
		double zeroVal = -9999999999;
		if (k == 0)
		{
			zeroVal = 0;
		}
		PureData[k][0] = zeroVal;


		for (int i = 1; i < Resolution; ++i)
		{
			double mu = i * dMu;

			double n = mu*mu/(var);
			double p = mu/(mu + var);
			double firstTerm = lgamma(n + k) - lgamma(n) - LogFac(k);
			double secondTerm = n * log(p) + k * log(1.0 - p);

			PureData[k][i] = (firstTerm + secondTerm);
		}
	}
}

