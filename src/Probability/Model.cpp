#include "Model.h"

void ModelParameters::Transform(const OptimiserParameters &in, int Kmax)
{
	Nu = 1+exp(in.x);
	Variance= 1+exp(in.y);
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

	Eta = 0.1/(1 + exp(-in.h));

	double Enorm = -9e99;
	int standardWidth = (Kmax+1)/in.gamma.size();
	double lWidth = log(standardWidth);
	int residualWidth = (Kmax + 1 - standardWidth * (in.gamma.size()-1));
	double rWidth = log(residualWidth);

	if (E.size() != in.gamma.size())
	{
		E.resize(in.gamma.size());
	}
	for (int i = 0; i < in.gamma.size(); ++i)
	{
		double width = lWidth;
		if (i == in.gamma.size() -1)	
		{
			width = rWidth;
		}
		E[i] = in.gamma[i];
		Enorm= ale(Enorm,E[i] + width);
	}

	for (int i = 0; i < in.gamma.size(); ++i)
	{
		E[i] -= Enorm;
	}

}

Model::Model(int kmax, int Q, int S, int errorRes)
{
	ErrorRes = errorRes;
	Kmax = kmax;
	NHarmonic = Q;
	Sum = S;
	ProbabilityArray.resize(kmax+1);
	Normalisation.resize(Q);
	logB.resize(Q,std::vector<double>(kmax+1,0.0));
	logK.resize(kmax+1);
	double prev = 0;
	logK[0] = 0;
	for (int k = 1; k <=kmax; ++k)
	{
		prev += log(k);
		logK[k] = prev;
	}

	SetParameters(OptimiserParameters(Q,errorRes));
	LOG(INFO) << "Probability Model Initialised with dimensions " << Q << "x" << kmax+1;
}

void Model::SetParameters(const OptimiserParameters & input)
{
	Parameters.Transform(input,Kmax);
	Compute();
}

double Model::LogError(int k)
{
	int rho = min((int)(k *1.0/(Kmax+1)* (ErrorRes)),ErrorRes-1);
	return Parameters.E[rho];
}
void Model::Compute()
{
	double snorm = -9999;
	for (int q = 0; q < NHarmonic; ++q)
	{
		double muq = (q + Parameters.Contamination[q]) * Parameters.Nu;
		double r = muq*muq/Parameters.Variance * Sum;
		double p = muq/(muq + Parameters.Variance);
		double rlogp = r*log(p);
		double log1mp = log1p(-p);
		double logW = log(Parameters.Weight[q]+1e-100);

		//accumulates logs as an efficient trick for computing lgamma(k+r) - lgamma(r)
		double logsum = 0;
		for (int k =0; k<= Kmax; ++k)
		{
			// double logNB = lgamma(k + r) - lgamma(r) - logK[k] + rlogp + k * log1mp;
			double logNB = logsum - logK[k] + rlogp + k*log1mp;
			logsum += log(k+r);

			logB[q][k] = logNB;
			if (k == 0)
			{
				Normalisation[q] = logNB;
			}
			else
			{
				Normalisation[q] = ale(Normalisation[q],logNB);
			}
			// double contribution = logW + logNB;
			// if (q == 0)
			// {
			// 	ProbabilityArray[k] = contribution;
			// }
			// else
			// {
			// 	ProbabilityArray[k] = ale(contribution,ProbabilityArray[k]);
			// }
			// snorm = ale(snorm,contribution);
		}
	}

	double sigFrac = log1p(-Parameters.Epsilon);
	double eFrac = log(Parameters.Epsilon);
	double s = 0;
	for (int k = 0; k<= Kmax; ++k)
	{
		double contribution = -9999;
		for (int q = 0; q <NHarmonic; ++q)
		{
			logB[q][k] -= Normalisation[q];
			contribution = ale(contribution,Parameters.LogWeight[q] + logB[q][k] );
		}
		// Normalisation = ale(Normalisation,ProbabilityArray[k]);
		ProbabilityArray[k] = ale(sigFrac + contribution, eFrac + LogError(k));
	}
}

double Model::Score(const std::vector<int> & histogram)
{
	double score = Prior();
	for (int k = 0; k < ProbabilityArray.size(); ++k)
	{
		score += histogram[k] * (ProbabilityArray[k]);
	}
	return score;
}

double Model::Prior()
{
	// return 0;
	double base = -1e2*pow(Parameters.Nu-83.0/Settings.Ploidy,2);
	// return base;
	double wplo = Parameters.Weight[Settings.Ploidy];
	for (int q = 0; q < Parameters.Weight.size(); ++q)
	{
		if (Parameters.Weight[q] > wplo)
		{
			base -= 1e5*pow(Parameters.Weight[q] - wplo,2);
		}

		base -= 1e3*pow(Parameters.Contamination[q] - (Settings.Ploidy - q)*Parameters.Eta,2);
		if (q > 0)
		{
			double sep = (Parameters.Contamination[q-1] - Parameters.Contamination[q]);
			double crit = 1.2;
			if (sep > crit)
			{
				base -= 1e5*pow(sep-crit,2);
			}
		}
	}
	return base;
}