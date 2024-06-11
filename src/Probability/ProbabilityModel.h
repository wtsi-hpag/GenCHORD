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
		ProbabilityModel(int NoiseResolution);
		ProbabilityModel(int NoiseResolution,int kMax, int Qmax);
		void SetDimensionality(int kMax, int Qmax);
		Dual GetDimensionality();
		void SetSignalParameters(double mu, double sigma);
		void SetSignalParameters(Dual input);


		
	
		double getLogNoise(int k);
		double getLogNoise(int k, int r);
		double getLogSignal(int k,int q);
		double logP(int k, int q);
		
		void SetGrids();


		void GlobalLogPrediction(std::vector<double> & output, const std::vector<double> & weights);

		Data Draw(int n,std::vector<double> ws, int kMax)
		{
			SetDimensionality(kMax,ws.size()-1);
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


		
		double SignalMean;
		double SignalSigma;

		double logNoiseWeight;
		double logSignalWeight;

		int NoiseResolution;
		std::vector<double> NoiseComponents;
		std::vector<double> logNoiseWidth;

		std::vector<double> Contamination;

		void SetNoiseComponents(std::vector<double> & wr);
		void ClearContamination();
		void SetContamination(int q,double targetMean);
		void ResetNoise();
		int RFromK(int k);
		// std::vector<double> KFromR(int r);
		int MaxK;
	protected:
		virtual double logNoise(int k);
		virtual double logSignal(int k, int q);
		int MaxQ;


		
		void SetNoiseGrid();
		void SetSignalGrids();


		std::vector<double> Noise;
		std::vector<std::vector<double>> Signal;

};





namespace Models
{
	// const double halflog2pi = 0.5 * log(2 * M_PI);
	// class Gaussian : public ProbabilityModel
	// {
	// 	public:
	// 		Gaussian() : ProbabilityModel(){};
	// 		Gaussian(double mu, double sigma,double eta, double mu_e,double sigma_e): ProbabilityModel()
	// 		{
	// 			SetSignalParameters(mu,sigma);
	// 			SetNoiseParameters(mu_e,sigma_e,exp(eta));
	// 		}

	// 		double logNoise(int k)
	// 		{
	// 			double mu_noise = NoiseMean;
	// 			double sigma_noise =NoiseSigma;
	// 			return -halflog2pi -log(sigma_noise)- 0.5 * pow((k - mu_noise)/sigma_noise,2);
	// 		}

	// 		double logSignal(int k, int q)
	// 		{
	// 			double mu = q * SignalMean;
	// 			double sigma = SignalSigma;
			
	// 			return -halflog2pi -log(sigma)- 0.5 * pow((k - mu)/sigma,2);

	// 		}

	// };

	class NegativeBinomial : public ProbabilityModel
	{
		public:
			NegativeBinomial(int NoiseResolution) : ProbabilityModel(NoiseResolution){};
			

			double logNoise(int k)
			{

				int r = RFromK(k);

				return NoiseComponents[r]- logNoiseWidth[r];
				// double mu = NoiseMean;
				// double sigma =NoiseSigma;
				// // std::cout << mu_noise << "  " << sigma_noise << std::endl;
				// double v = sigma*sigma;
				// double r = mu*mu/v;
				// double p = mu/(mu + v);

				// double t1 = k * log(1 - p) + r * log(p);
				// double t2 = lgamma(k+r) - lgamma(k+1) - lgamma(r);

				// return t1+t2;
			}

			double logSignal(int k, int q)
			{
				double c= Contamination[q];
				double mu = (q * (1.0 - c) + 2 *c)* SignalMean;
				mu = std::max(0.1,mu);
				double sigma = SignalSigma;
				// if (q == 0)
				// {
				// 	mu = 1;
				// 	// mu = std::max(2);//std::min(sigma/6,SignalMean/8);
				// 	// sigma /= 2;
				// }
				// if (q == 1)
				// {
				// 	mu*=1.2;
				// }
				double v = sigma*sigma;
				double r = mu*mu/v;
				double p = mu/(mu + v);

				double t1 = k * log(1 - p) + r * log(p);
				double t2 = lgamma(k+r) - lgamma(k+1) - lgamma(r);

				return t1+t2;
			}

	};
}
