#pragma once
#include "../Utility/basicFunctions.h"
#include "../DataFrame/data.h"
#include <vector>
#include <iostream>

struct Dual
{
	double X;
	double Y;
	Dual(){
		X = 0;
		Y = 0;
	}
	Dual(double x, double y)
	{
		X = x;
		Y = y;
	}
};

class ProbabilityModel
{
	public:
		ProbabilityModel();
		ProbabilityModel(int kMax, int Qmax);
		void SetDimensionality(int kMax, int Qmax);
		Dual GetDimensionality();
		void SetSignalParameters(double mu, double sigma);
		void SetSignalParameters(Dual input);
		void SetNoiseParameters(double mu, double sigma, double gamma);
		void SetNoiseParametersFromObserved(double mu,double sigma, double gamma);


		
	
		double getLogNoise(int k);
		double getLogSignal(int k,int q);
		double logP(int k, int q);
		
		void SetGrids();


		void GlobalLogPrediction(std::vector<double> & output, const std::vector<double> & weights);

		Data Draw(int n,std::vector<double> ws, int kMax)
		{
			SetDimensionality(kMax,ws.size()-1);
			std::cout << MaxK << "  " << MaxQ << std::endl;
			ws = JSL::Vector(ws)/JSL::Vector(ws).Sum();
			std::vector<int> ks = JSL::Vector::intspace(0,kMax,1);

			SetGrids();
			std::vector<double> probs(ks.size());
			GlobalLogPrediction(probs,ws);

			std::vector<double> cumsum(probs.size(),0);
			cumsum[0] = exp(probs[0]);
			for (int i = 1; i < ks.size(); ++i)
			{
				cumsum[i] = cumsum[i-1] + exp(probs[i]);
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


		int NoiseResolution = 180;
		double NoiseMeanMin=1;
		double NoiseMeanMax = 150;
		double NoiseSigmaMin = 0.01;
		double NoiseSigmaMax = 70;

		double SignalMean;
		double SignalSigma;

		double NoiseWeight;
		double NoiseMean;
		double NoiseSigma;
		double logNoiseWeight;
		double logSignalWeight;

		void ComputeNoiseConversionChart();
	protected:
		virtual double logNoise(int k);
		virtual double logSignal(int k, int q);
		int MaxK;
		int MaxQ;

		
		double inferredMeanMin;
		double inferredMeanMax;
		double inferredSigmaMin;
		double inferredSigmaMax;
		
		
		void SetNoiseGrid();
		void SetSignalGrids();
		std::vector<std::vector<Dual>> NoiseGrid;

		Dual GetMuSigma(double mean, double sigma);

		Dual ComputeNoiseMoments();


		std::vector<double> Noise;
		std::vector<std::vector<double>> Signal;

};





namespace Models
{
	const double halflog2pi = 0.5 * log(2 * M_PI);
	class Gaussian : public ProbabilityModel
	{
		public:
			Gaussian() : ProbabilityModel(){};
			Gaussian(double mu, double sigma,double eta, double mu_e,double sigma_e): ProbabilityModel()
			{
				SetSignalParameters(mu,sigma);
				SetNoiseParameters(mu_e,sigma_e,exp(eta));
			}

			double logNoise(int k)
			{
				double mu_noise = NoiseMean;
				double sigma_noise =NoiseSigma;
				return -halflog2pi -log(sigma_noise)- 0.5 * pow((k - mu_noise)/sigma_noise,2);
			}

			double logSignal(int k, int q)
			{
				double mu = q * SignalMean;
				double sigma = SignalSigma;
			
				return -halflog2pi -log(sigma)- 0.5 * pow((k - mu)/sigma,2);

			}

	};

	class NegativeBinomial : public ProbabilityModel
	{
		public:
			NegativeBinomial() : ProbabilityModel(){};
			NegativeBinomial(double mu, double sigma,double eta, double mu_e,double sigma_e): ProbabilityModel()
			{
				SetSignalParameters(mu,sigma);
				SetNoiseParameters(mu_e,sigma_e,exp(eta));
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
					mu = 1;
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
