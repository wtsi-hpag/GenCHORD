#include "Model.h"

void ModelParameters::Transform(const OptimiserParameters &in)
{
	Nu = exp(in.x);
	Variance= exp(in.y);
	Epsilon = Settings.ErrorMax/ ( 1 + exp(-in.phi));
	
	double s = 0;
	int Q = Weight.size();
	for (int q= 0; q < Q; ++q)
	{
		auto expz = exp(in.z[q]);
		s += expz;
		Weight[q] = s;

		double dminus = 0;
		if (q > 0)
		{
			dminus = Settings.ContaminationMin;
		}
		double dplus = Settings.ContaminationMax;
		Contamination[q] = dminus + (dplus - dminus)/(1 + exp(-in.psi[q]));
	}

	for (int q= 0; q < Q; ++q)
	{
		Weight[q]/= s;

	}
}

Model::Model(int kmax, int Q, int S)
{
	Kmax = kmax;
	NHarmonic = Q;
	Sum = S;
	ProbabilityArray.resize(kmax+1);
	logK.resize(kmax+1);
	double prev = 0;
	logK[0] = 0;
	for (int k = 1; k <=kmax; ++k)
	{
		prev += log(k);
		logK[k] = prev;
	}
}

double Model::LogError(int k)
{
	return -log(Kmax);
}
void Model::Compute()
{
	
	for (int q = 0; q < NHarmonic; ++q)
	{
		double muq = (q + Parameters.Contamination[q]) * Parameters.Nu;
		double r = muq*muq/Parameters.Variance * Sum;
		double p = muq/(muq + Parameters.Variance);
		LOG(DEBUG) << r << " " << p;
		double logp = log(p);
		double log1mp = log1p(-p);
		double logW = log(Parameters.Weight[q]+1e-100);
		for (int k =0; k<= Kmax; ++k)
		{
			double logNB = lgamma(k + r) - lgamma(r) - logK[k] + r * logp + k * log1mp;

			if (q == 0)
			{
				ProbabilityArray[k] = logW + logNB;
			}
			else
			{
				ProbabilityArray[k] = ale(logW + logNB,ProbabilityArray[k]);
			}
		}
	}

	double sigFrac = log1p(-Parameters.Epsilon);
	double eFrac = log(Parameters.Epsilon);
	double s = 0;
	for (int k = 0; k<= Kmax; ++k)
	{
		ProbabilityArray[k] = ale(sigFrac + ProbabilityArray[k], eFrac + LogError(k));
		if (k == 0)
		{
			s = ProbabilityArray[k];
		}
		else
		{
			s = ale(ProbabilityArray[k],s);
		}
	}
	for (int k = 0; k<= Kmax; ++k)
	{
		ProbabilityArray[k] -= s;
	}
}