#pragma once
#include "../settings.h"
#include "../Utility/ale.h"
#include "../Utility/Random.h"
static const double v = 1.0/RAND_MAX;

struct OptimiserParameters
{
	double x;
	double y;
	std::vector<double> z;
	double phi;
	std::vector<double> psi;
	std::vector<double> gamma;
	double h;
	OptimiserParameters(int dim=0, int res = 20)
	{
		x = log(30);
		y = 5;
		z = std::vector<double>(dim,0.0);
		psi = std::vector<double>(dim,0.01);
		phi =0;
		h = 0;
		gamma = std::vector<double>(res,0);
		// LOG(ERROR)
	}
	template <typename T>
	void RandomStep(Random<T> & R, OptimiserParameters & propose,double s = 1)
	{
		// double s = 1;
		propose.x = x + R.Normal(0,s);
		propose.y = y + R.Normal(0,s);
		propose.phi = phi + R.Normal(0,s);
		propose.h =  h + R.Normal(0,s);
		for (int i = 0; i < z.size(); ++i)
		{
			propose.z[i] = z[i] + R.Normal(0,s);
			propose.psi[i] = psi[i] + R.Normal(0,s);
		}
		for (int rho = 0; rho < gamma.size(); ++rho)
		{
			propose.gamma[rho] = gamma[rho] + R.Normal(0,s);
		}
		propose.psi[Settings.Ploidy] = -100;

	}
};

struct ModelParameters
{
	double Nu;
	double Variance;
	std::vector<double> Weight;
	std::vector<double> LogWeight;
	double Epsilon;
	std::vector<double> Contamination;
	double Eta;
	std::vector<double> E;
	void Transform(const OptimiserParameters & in, int kmax);
};

class Model
{
	private:
		std::vector<double> ProbabilityArray;
		std::vector<double> logK;
		std::vector<double> Normalisation;
		
		std::vector<std::vector<double>>logB;
		double eNorm;
		int ErrorRes;
	public: 
		int Kmax;
		int NHarmonic;
		int Sum;
		ModelParameters Parameters;
		Model(int kmax, int Q, int S, int res);
	
		double Prior();
		double Score(const std::vector<int> & histogram);
		double LogError(int k);
		void Compute();
		double operator[](int s) const {return ProbabilityArray[s];}; // deliberately does not have the non-const version. Altering the probability array can only be done by updating the parameters.
		void SetParameters(const OptimiserParameters & input);
		double Sample(int q, int s)  {return Parameters.LogWeight[q] + logB[q][s] + log1p(-Parameters.Epsilon);}
};