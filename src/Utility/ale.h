#pragma once


//computes log(e^x + e^y) in an overflow-safe fashion
//use of log1p increases precision/performance.
inline double ale(double x, double y)
{
	return max(x,y) + log1p(exp(-abs(x-y)));
}