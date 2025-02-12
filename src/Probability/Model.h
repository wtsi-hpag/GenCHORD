#pragma once
#include "../settings.h"
#include "../Utility/ale.h"
struct OptimiserParameters
{
	double x;
	double y;
	std::vector<double> z;
	double phi;
	std::vector<double> psi;
};

struct ModelParameters
{
	double Nu;
	double Variance;
	std::vector<double> Weight;
	double Epsilon;
	std::vector<double> Contamination;
	void Transform(const OptimiserParameters & in);
};

class Model
{
	private:
		std::vector<double> ProbabilityArray;
		std::vector<double> logK;
	public: 
		int Kmax;
		int NHarmonic;
		int Sum;
		ModelParameters Parameters;
		Model(int kmax, int Q, int S);
	
		double LogError(int k);
		void Compute();
		double operator[](int s) const {return ProbabilityArray[s];};
		double & operator[](int s) {return ProbabilityArray[s];};
};