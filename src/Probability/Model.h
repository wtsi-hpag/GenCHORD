#pragma once
#include "../settings.h"
#include "../Utility/ale.h"
#include "../Utility/Random.h"
#include "Parameters.h"


class Model
{
	private:
		std::vector<double> ProbabilityArray;
		std::vector<std::vector<double>> HarmonicProbabilityArray;
		std::vector<double> logK;
		
		std::vector<std::vector<double>>logB;
		double eNorm;
		int ErrorRes;
		int ErrorWindow(int k);
		std::vector<double> sMean;
		std::vector<double> digammaMean;
		std::vector<double> digammaArray;
		std::vector<double> contaminationPref;

		std::vector<double> Muq;
		std::vector<double> Rq;
		std::vector<double> Pq;
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
		void SetParameters(const OptimiserPack & input);
		void SetParameters(const StateVector & input);
		void ComputeGradient(StateVector & grad,const std::vector<int> & histogram);
		double Sample(int q, int s)  {return Parameters.LogWeight[q] + logB[q][s] + log1p(-Parameters.Epsilon);}
		void PrepareHarmonics();
		double HarmonicProbability(int k, int q);
};