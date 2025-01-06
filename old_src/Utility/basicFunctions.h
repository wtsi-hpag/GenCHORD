#pragma once
#include <math.h>
// #include <algorithm>
#include <vector>
extern bool globalVerbose;
#define Log( ...)			\
{									\
	if (globalVerbose) 					\
	{								\
		std::cout << __VA_ARGS__;	\
	}								\
}

double ale(double x, double y);


std::vector<double> exp(std::vector<double> x);