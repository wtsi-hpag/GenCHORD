#pragma once
#include "../Utility/basicFunctions.h"
#include "../DataFrame/data.h"
#include <vector>
#include <iostream>

struct MuSigmaPair
{
	double Mu;
	double Sigma;
	MuSigmaPair(){
		Mu = 0;
		Sigma = 0;
	}
	MuSigmaPair(double m, double s)
	{
		Mu = m;
		Sigma = s;
	}
};

class ProbabilityModel
{
	public:
		int Dimensionality;
		ProbabilityModel(int dimensions);

		virtual double logNoise(int k);

		virtual double logSignal(int k, int q);
	
		double logP(int k, int q,bool withCorrections);

		double logP(int k, int q);
		void FillGrid(std::vector<std::vector<double>> & grid);

		void SetParameters(std::vector<double> v);

		void Normalise(int kMax, std::vector<double> & weights);

		std::vector<double> GlobalPrediction(int kMax, std::vector<double> & weights);
		std::vector<double> Parameters;

		Data Draw(int n,std::vector<double> ws, int kMax)
		{
			ws = JSL::Vector(ws)/JSL::Vector(ws).Sum();
			std::vector<int> ks = JSL::Vector::intspace(0,kMax,1);

			auto probs = GlobalPrediction(kMax,ws);
			std::vector<double> cumsum(probs.size(),0);
			cumsum[0] = probs[0];
			for (int i = 1; i < ks.size(); ++i)
			{
				cumsum[i] = cumsum[i-1] + probs[i];
			}
			double fin = cumsum[ks.size()-1];
			ChromosomeCoverage c;

			double s = 0;

			for (int i = 0; i < n; ++i)
			{
				double r = rand()*1.0/RAND_MAX;
				int chosen = -1;
				for (int j = 0; j < cumsum.size();++j)
				{
					if (r < cumsum[j]/fin)
					{
						c.Add(i,ks[j]);
						s+= ks[j];
						break;
					}
				}
			}

			Data d;
			d.Mean = s/n;
			d.Chromosomes = {c};
			return d;
		}



		double NoiseNormalisation;
		std::vector<double> Normalisation;

		int NoiseResolution = 150;
		double NoiseMeanMin=1;
		double NoiseMeanMax = 170;
		double NoiseSigmaMin = 0.01;
		double NoiseSigmaMax = 60;

		double SignalMean;
		double SignalSigma;

		double NoiseWeight;
		double NoiseMean;
		double NoiseSigma;

		void SetNoiseParameters_Observed(double weight, double observedMean, double observedSigma);
		void ComputeNoiseConversionChart(int kMax);
	protected:

		double inferredMeanMin;
		double inferredMeanMax;
		double inferredSigmaMin;
		double inferredSigmaMax;
		std::vector<std::vector<MuSigmaPair>> NoiseGrid;

		MuSigmaPair GetMuSigma(double mean, double sigma);

		MuSigmaPair ComputeMoments(int kMax);

};





namespace Models
{
	const double halflog2pi = 0.5 * log(2 * M_PI);
	class Gaussian : public ProbabilityModel
	{
		public:
			Gaussian() : ProbabilityModel(2){};
			Gaussian(double mu, double sigma,double eta, double mu_e,double sigma_e): ProbabilityModel(5)
			{
				SignalMean = mu;
				SignalSigma = sigma;
				NoiseWeight = exp(eta);
				NoiseMean = mu_e;
				NoiseSigma = sigma_e;
			}

			double logNoise(int k)
			{
				double mu_noise = NoiseMean;
				double sigma_noise =NoiseSigma;
				// std::cout << mu_noise << "  " << sigma_noise << std::endl;
				return -halflog2pi -log(sigma_noise)- 0.5 * pow((k - mu_noise)/sigma_noise,2);
			}

			double logSignal(int k, int q)
			{
				double mu = q * SignalMean;
				double sigma = SignalSigma;
				if (q == 0)
				{
					sigma/=2;
				}
				return -halflog2pi -log(sigma)- 0.5 * pow((k - mu)/sigma,2);

			}

	};

	class NegativeBinomial : public ProbabilityModel
	{
		public:
			NegativeBinomial() : ProbabilityModel(2){};
			NegativeBinomial(double mu, double sigma,double eta, double mu_e,double sigma_e): ProbabilityModel(5)
			{
				SignalMean = mu;
				SignalSigma = sigma;
				NoiseWeight = exp(eta);
				NoiseMean = mu_e;
				NoiseSigma = sigma_e;
			}

			double logNoise(int k)
			{
				double mu = NoiseMean;
				double sigma =NoiseSigma;
				// std::cout << mu_noise << "  " << sigma_noise << std::endl;
				double v = sigma*sigma;
				double r = mu*mu/v;
				double p = mu/(mu + v);

				double t1 = k * log(1 - p) + r * log(p);
				double t2 = lgamma(k+r) - lgamma(k+1) - lgamma(r);

				return t1+t2;
			}

			double logSignal(int k, int q)
			{
				double mu = q * SignalMean;
				double sigma = SignalSigma;
				if (q == 0)
				{
					mu = sigma/10;
				}
				double v = sigma*sigma;
				double r = mu*mu/v;
				double p = mu/(mu + v);

				double t1 = k * log(1 - p) + r * log(p);
				double t2 = lgamma(k+r) - lgamma(k+1) - lgamma(r);

				return t1+t2;
			}

	};
}
