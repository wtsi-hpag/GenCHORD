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

	OptimiserParameters(int dim)
	{
		x = log(30);
		y = 0;
		z = std::vector<double>(dim,0.0);
		psi = std::vector<double>(dim,0.01);
	}
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

struct lgammaLookup
{
	int Resolution;
	int Maximum;
	double Delta;
	std::vector<double> values;
	double operator[](double s) const;
	void Initialise(int max,int res);
};

class Model
{
	private:
		std::vector<double> ProbabilityArray;
		std::vector<double> logK;
		std::vector<std::vector<double>>logB;
		lgammaLookup logGamma;
	public: 
		int Kmax;
		int NHarmonic;
		int Sum;
		ModelParameters Parameters;
		Model(int kmax, int Q, int S);
	
		double Prior();
		double Score(const std::vector<int> & histogram);
		double LogError(int k);
		void Compute();
		double operator[](int s) const {return ProbabilityArray[s];};
		double & operator[](int s) {return ProbabilityArray[s];};
		void SetParameters(const OptimiserParameters & input);
};