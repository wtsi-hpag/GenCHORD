#pragma once
#include <vector>
#include "../Settings.h"
#include "../Utility/Random.h"
class StateVector
{
	private:
		int AccumulateOrder;
	public:
		double x;
		double y;
		std::vector<double> z;
		double phi;
		std::vector<double> psi;
		std::vector<double> gamma;
		double h;

		StateVector(int dim, int res, int accumulateOrder=0);

		//templated to allow different random engines to be used.
		template <typename T>
		void RandomStep(Random<T> & R, StateVector & propose,double step = 1)
		{
			propose.x = x + R.Normal(0,step);
			propose.y = y + R.Normal(0,step);
			propose.phi = phi + R.Normal(0,step);
			propose.h =  h + R.Normal(0,step);
			for (int i = 0; i < z.size(); ++i)
			{
				propose.z[i] = z[i] + R.Normal(0,step);
				propose.psi[i] = psi[i] + R.Normal(0,step);
			}
			for (int rho = 0; rho < gamma.size(); ++rho)
			{
				propose.gamma[rho] = gamma[rho] + R.Normal(0,step);
			}
			propose.psi[Settings.Ploidy] = -100;

		}

		void SetDefaultValues();
		void Accumulate(const StateVector & gradient, double memory);

		void ADAMUpdate(const StateVector & firstMoment, const StateVector & secondMoment, double alpha, double b1, double b2, int l);

		void Reset();

		void LightCopy(const StateVector & in);
};