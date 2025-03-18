#include "basicFunctions.h"

bool globalVerbose = true;

double ale(double x, double y)
{
	return std::max(x,y) + log(1.0 + exp(-abs(x - y)));
}

std::vector<double> exp(std::vector<double> x)
{
	std::vector<double> out(x.size());
	for (int i =0; i < x.size(); ++i)
	{
		out[i] = exp(x[i]);
	}
	return out;
}