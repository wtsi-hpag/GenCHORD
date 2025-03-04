#include "Model.h"


Model::Model(int kmax, int Q, int S, int errorRes)
{
	ErrorRes = errorRes;
	Kmax = kmax;
	NHarmonic = Q;
	Sum = S;
	ProbabilityArray.resize(kmax+1);
	logB.resize(Q,std::vector<double>(kmax+1,0.0));
	logK.resize(kmax+1);
	double prev = 0;
	logK[0] = 0;
	for (int k = 1; k <=kmax; ++k)
	{
		prev += log(k);
		logK[k] = prev;
	}

	SetParameters(OptimiserPack(Q,errorRes));
	LOG(INFO) << "Probability Model Initialised with dimensions " << Q << "x" << kmax+1;
}

void Model::SetParameters(const StateVector & input)
{
	Parameters.Transform(input,Kmax);
	Compute();
}
void Model::SetParameters(const OptimiserPack & input)
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

	//compute B_k^q, the negative binomial grid
	for (int q = 0; q < NHarmonic; ++q)
	{
		//NB parameters
		double muq = (q + Parameters.Contamination[q]) * Parameters.Nu;
		double r = muq*muq/Parameters.Variance * Sum;
		double p = muq/(muq + Parameters.Variance);
		
		//some shortcuts to avoid recalculating at every step
		double rlogp = r*log(p);
		double log1mp = log1p(-p);
		double logW = log(Parameters.Weight[q]+1e-100);

		//accumulates logs as an efficient trick for computing lgamma(k+r) - lgamma(r)
		double logsum = 0;

		//require logNB to sum to one across domain
		double normSum = -999999;

		for (int k =0; k<= Kmax; ++k)
		{
			double logNB = logsum - logK[k] + rlogp + k*log1mp;
			logsum += log(k+r); //<- this is the trick discussed above

			logB[q][k] = logNB;
			normSum = ale(normSum,logNB);
		}

		//normalise B_k^q to one across k.
		for (int k =0; k<= Kmax; ++k)
		{
			logB[q][k] -= normSum;
		}
	}

	//add in the error component.
	double sigFrac = log1p(-Parameters.Epsilon);
	double eFrac = log(Parameters.Epsilon);
	double s = 0;
	for (int k = 0; k<= Kmax; ++k)
	{
		double contribution = eFrac + LogError(k);
		for (int q = 0; q <NHarmonic; ++q)
		{
			contribution = ale(contribution,Parameters.LogWeight[q] + logB[q][k] );
		}
		ProbabilityArray[k] = contribution;
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
	// double base = 0;
	// if (Parameters.Nu > Kmax * 1.0/Settings.AccumulationFactor)
	// {
	// 	base -= pow(base - Kmax *1.0/Settings.AccumulationFactor,2);
	// }
	// return base;
	double base = -1e1*pow(Parameters.Nu-83.0/Settings.Ploidy,2);
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