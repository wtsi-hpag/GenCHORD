#include "Parameters.h"

OptimiserPack::OptimiserPack(int harmonicCount, int errorResolution): Parameters(harmonicCount,errorResolution), Gradient(harmonicCount, errorResolution), FirstMoment(harmonicCount,errorResolution,1), SecondMoment(harmonicCount,errorResolution,2)
{
	Parameters.SetDefaultValues();
}

OptimiserPack::OptimiserPack(StateVector input): Parameters(input), Gradient(input.z.size(),input.gamma.size()), FirstMoment(input.z.size(),input.gamma.size(),1), SecondMoment(input.z.size(),input.gamma.size(),2)
{

}



void OptimiserPack::AccumulateGradient(double b1, double b2)
{
	FirstMoment.Accumulate(Gradient,b1);
	SecondMoment.Accumulate(Gradient,b2);
}
void OptimiserPack::ADAMUpdate(double alpha, double b1, double b2, int l)
{
	Parameters.ADAMUpdate(FirstMoment,SecondMoment,alpha,b1,b2,l);
}

void OptimiserPack::OptimiseReset()
{
	Gradient.Reset();
	FirstMoment.Reset();
	SecondMoment.Reset();
}

void ModelParameters::Transform(const OptimiserPack &in, int Kmax)
{
	Transform(in.Parameters, Kmax);
}

void ModelParameters::SetWindows(int Kmax)
{
	StandardWidth = floor((Kmax+1)/E.size() + 0.5);
	logWidth = log(StandardWidth);
	ResidualWidth = (Kmax + 1 - StandardWidth * (E.size()-1));
	logRWidth = log(ResidualWidth);
	
	WindowSet = true;
}

void ModelParameters::Transform(const StateVector &in, int Kmax)
{
	Nu = exp(in.x);
	Variance= exp(in.y);
	Epsilon = Settings.ErrorMax/ ( 1 + exp(-in.phi));
	double s = 0;
	int Q = in.z.size();
	if (Q != Weight.size())
	{
		Weight.resize(Q,0.0);
		Contamination.resize(Q,0.0);
		LogWeight.resize(Q,0.0);
	}
	for (int q= 0; q < Q; ++q)
	{
		auto expz = exp(in.z[q]);
		s += expz;
		Weight[q] = expz;
		
		double dminus = 0;
		if (q > 0)
		{
			dminus = Settings.ContaminationMin;
		}
		double dplus = Settings.ContaminationMax;
		Contamination[q] = dminus + (dplus - dminus)/(1 + exp(-in.psi[q]));
	}
	Contamination[Settings.Ploidy] =0;

	for (int q= 0; q < Q; ++q)
	{
		Weight[q]/= s;
		LogWeight[q] = log(Weight[q]);
	}

	Eta = Settings.ContaminationMax/(1 + exp(-in.h));

	double Enorm = -9e99;

	if (E.size() != in.gamma.size() || !WindowSet)
	{
		E.resize(in.gamma.size());
		SetWindows(Kmax);
	}
	for (int i = 0; i < in.gamma.size(); ++i)
	{
		double width = logWidth;
		if (i == in.gamma.size() -1)	
		{
			width = logRWidth;
		}
		E[i] = in.gamma[i];
		Enorm= ale(Enorm,E[i] + width);
		
	}

	for (int i = 0; i < in.gamma.size(); ++i)
	{
		E[i] -= Enorm;
	}
}