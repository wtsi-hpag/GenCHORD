#include "StateVector.h"
StateVector::StateVector(int dim, int res, int acc)
{
	//this initialises to the `zero vector' which is suitable for the default containers. The actual optimisation vector itself probably doesn't want to be zero initialised -- we provide a default value function to call.
	x = 0;
	y = 0;
	z = std::vector<double>(dim,0.0);
	psi = std::vector<double>(dim,0.0);
	if (dim <= Settings.Ploidy)
	{
		throw std::runtime_error("Cannot construct parameter vectors with Q value less than or equal to the Ploidy! You must include higher order harmonics for the analysis to function.");
	}
	phi =0;
	h = 0;
	gamma = std::vector<double>(res,0);
	AccumulateOrder = acc;
}

void StateVector::LightCopy(const StateVector & in)
{
	x = in.x;
	y = in.y;
	z = in.z;
	phi = in.phi;
	psi = in.psi;
	gamma = in.gamma;
	h = in.h;
}

void StateVector::SetDefaultValues()
{
	x = log(Settings.DefaultMean/Settings.Ploidy);
	y = 1;
	phi = -1;
	h = -1;
	std::fill(z.begin(), z.end(),0.0);
	std::fill(psi.begin(), psi.end(),0);
	std::fill(gamma.begin(),gamma.end(),0.0);
}

void StateVector::Accumulate(const StateVector & gradient, double memory)
{
	// int order = AccumulateOrder;
	double inv = 1.0 - memory;
	x = memory * x + inv * pow(gradient.x,AccumulateOrder);
	y = memory * y + inv * pow(gradient.y,AccumulateOrder);
	phi = memory * phi + inv * pow(gradient.phi,AccumulateOrder);
	h = memory * h + inv * pow(gradient.h,AccumulateOrder);

	for (int i = 0; i < z.size(); ++i)
	{
		z[i] = memory * z[i] + inv * pow(gradient.z[i],AccumulateOrder);
		psi[i] = memory * psi[i] + inv * pow(gradient.psi[i],AccumulateOrder);
	}
	for (int j = 0; j < gamma.size(); ++j)
	{
		gamma[j] = memory * gamma[j] + inv * pow(gradient.gamma[j],AccumulateOrder);
	}
}
void StateVector::ADAMUpdate(const StateVector & firstMoment, const StateVector & secondMoment, double alpha, double b1, double b2, int l)
{
	const double c1 = 1.0/(1.0 - pow(b1,l));
	const double c2 = 1.0/(1.0 - pow(b2,l));
	const double eps = 1e-10; 
	x += 0.1*alpha * firstMoment.x * c1 / sqrt(secondMoment.x*c2 + eps);
	y += alpha * firstMoment.y * c1 / sqrt(secondMoment.y*c2 + eps);
	phi += alpha * firstMoment.phi * c1 / sqrt(secondMoment.phi*c2 + eps);
	h += alpha * firstMoment.h * c1 / sqrt(secondMoment.h*c2 + eps);	
	
	for (int i = 0; i < z.size(); ++i)
	{
		z[i] += alpha * firstMoment.z[i] * c1/sqrt(secondMoment.z[i] + eps);
		psi[i] += alpha * firstMoment.psi[i] * c1 / sqrt(secondMoment.psi[i]*c2 + eps);
	}
	for (int j = 0; j < gamma.size(); ++j)
	{
		gamma[j] += alpha * firstMoment.gamma[j] * c1/sqrt(secondMoment.gamma[j]*c2 + eps);
	}

	if (x> log(150))
	{
		x = log(150);
	}
}

void StateVector::Reset()
{
	x = 0;
	y = 0;
	phi = 0;
	h = 0;
	std::fill(z.begin(), z.end(),0.0);
	std::fill(psi.begin(), psi.end(),0.0);
	std::fill(gamma.begin(),gamma.end(),0.0);
}
