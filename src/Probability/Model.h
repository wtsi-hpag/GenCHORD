#pragma once
#include "../settings.h"
#include "../Utility/ale.h"

static const double v = 1.0/RAND_MAX;
double inline random(double min, double max)
{
	return min + (max - min) * rand()*v;
}

struct OptimiserParameters
{
	double x;
	double y;
	std::vector<double> z;
	double phi;
	std::vector<double> psi;
	double h;
	OptimiserParameters(int dim=0)
	{
		x = log(30);
		y = 0;
		z = std::vector<double>(dim,0.0);
		psi = std::vector<double>(dim,0.01);
		phi =0;
		h = 0;
	}

	void Randomise()
	{
		x = random(log(10),log(35));
		y = random(0,8);
		phi = random(-5,5);
		h = random(-5,5);
		for (int i = 0; i < z.size(); ++i)
		{
			z[i] = random(-10,10);
			psi[i] = random(-10,10);
			psi[Settings.Ploidy] = -100;
		}
	}
};

struct ModelParameters
{
	double Nu;
	double Variance;
	std::vector<double> Weight;
	double Epsilon;
	std::vector<double> Contamination;
	double Eta;
	void Transform(const OptimiserParameters & in);
};

class Model
{
	private:
		std::vector<double> ProbabilityArray;
		std::vector<double> logK;
		std::vector<std::vector<double>>logB;
	public: 
		int Kmax;
		int NHarmonic;
		int Sum;
		double Normalisation;
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