#include "TruncatedGaussian.h"

TruncatedGaussian::TruncatedGaussian(int kMax, int bounder)
{
	Truncator = bounder;
	KMax = kMax;
	PureData.resize(kMax+1,std::vector<double>(kMax + bounder+2,0.0));
}
void TruncatedGaussian::Populate(double gamma)
{
	for (int kobs = 0; kobs <= KMax; ++kobs)
	{
		double s = -999999999;
		for (int k = std::max(0,kobs-Truncator); k <= KMax + Truncator; ++k)
		{
			double d = (k - kobs)/gamma;
			double v = - 0.5 * d * d;
			PureData[kobs][k] = v;
			s = ale(s,v);
		}
		for (int k = std::max(0,kobs-Truncator); k <= KMax + Truncator; ++k)
		{
			PureData[kobs][k] -= s;
		}
	}
}