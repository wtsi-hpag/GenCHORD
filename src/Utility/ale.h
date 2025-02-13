#pragma once



inline double tailoredFastLog1p(double x)
{
    if (x < -37.0) return exp(x);   // log(1+exp(x)) ≈ exp(x) for very small x. exp(-37) is machine epsilon, so 1 + 
    return log1p(exp(x));  	
}

//computes log(e^x + e^y) in an overflow-safe fashion
//use of log1p increases precision/performance.
inline double ale(double x, double y)
{
	return max(x,y) + tailoredFastLog1p(-abs(x-y));
	
	// log1p(exp(-abs(x-y)));
}