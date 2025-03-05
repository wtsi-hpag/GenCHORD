#pragma once
#include "../settings.h"
#include "../Utility/Random.h"
#include "../Utility/ale.h"
#include "StateVector.h"
class OptimiserPack
{
	public:
		
		StateVector Parameters;
		StateVector Gradient;

		OptimiserPack(StateVector input);
		OptimiserPack(int harmonicCount, int errorResolution);

		void AccumulateGradient(double b1, double b2);
		void ADAMUpdate(double alpha, double b1, double b2, int l);
		// void Accumulate(const OptimiserParameters & gradient, double memory, int order);
		void OptimiseReset();
	private:
		StateVector FirstMoment;
		StateVector SecondMoment;
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
	void SetWindows(int Kmax);
	void Transform(const OptimiserPack & in, int kmax);
	void Transform(const StateVector & in, int kmax);

	// private:
	int StandardWidth;
	double logWidth;
	int ResidualWidth;
	double logRWidth;
	bool WindowSet = false;
};