#include "Model.h"


double digammaApprox(double x)
{
	if (x < 5)
	{
		return digammaApprox(x+1) - 1.0/x;
	}
	return log(x) - 1./(2*x) - 1.0/(12 * x*x) + 1.0/(120 * pow(x,4));
}

Model::Model(int kmax, int Q, int S, int errorRes)
{
	ErrorRes = errorRes;
	Kmax = kmax;
	NHarmonic = Q;
	Sum = S;
	ProbabilityArray.resize(kmax+1);
	logB.resize(Q,std::vector<double>(kmax+1,0.0));
	digammaMean.resize(Q);
	digammaArray.resize(Q);
	contaminationPref.resize(Q);
	Muq.resize(Q);
	Rq.resize(Q);
	Pq.resize(Q);
	sMean.resize(Q);
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

int Model::ErrorWindow(int k)
{
	if (k < Kmax +1 - Parameters.StandardWidth)
	{
		return k/Parameters.StandardWidth;
	}
	else
	{
		return ErrorRes - 1;
	}
}

double Model::LogError(int k)
{
	return Parameters.E[ErrorWindow(k)];
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
		Muq[q] = muq;
		Rq[q] = r;
		Pq[q] = p;
		
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
		sMean[q] = 0;
		digammaMean[q] = 0;
		double digamma = 0;
		for (int k =0; k<= Kmax; ++k)
		{
			logB[q][k] -= normSum;
			double temp = exp(logB[q][k]);
			sMean[q] += k*temp;
			digammaMean[q] += digamma * temp; 
			digamma += 1.0/(k+r);
		}
		// LOG(DEBUG) << "\tqMean_" << q << " = " << sMean[q] << " " << digammaMean[q] << " (from mu = " << muq <<")" << " " << normSum;
	}

	

	//add in the error component.
	double sigFrac = log1p(-Parameters.Epsilon);
	double eFrac = log(Parameters.Epsilon);
	
	for (int k = 0; k<= Kmax; ++k)
	{
		double contribution = eFrac + LogError(k);
		for (int q = 0; q <NHarmonic; ++q)
		{
			contribution = ale(contribution,sigFrac + Parameters.LogWeight[q] + logB[q][k] );
		}
		ProbabilityArray[k] = contribution;
	}
}

void Model::PrepareHarmonics()
{
	HarmonicProbabilityArray.resize(Kmax+1,std::vector<double>(NHarmonic,0.0));

	double sigFrac = log1p(-Parameters.Epsilon);
	double eFrac = log(Parameters.Epsilon);
	
	for (int k = 0; k<= Kmax; ++k)
	{
		double error = eFrac + LogError(k);
		for (int q = 0; q <NHarmonic; ++q)
		{
			HarmonicProbabilityArray[k][q] = ale(error,sigFrac + logB[q][k]);
		}
	}
}
double Model::HarmonicProbability(int k, int q)
{
	if (k > Kmax)
	{
		k = Kmax;
	}
	return HarmonicProbabilityArray[k][q];
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
	// double base = 0;
	// if (Parameters.Nu > Kmax * 1.0/Settings.AccumulationFactor)
	// {
	// 	base -= pow(base - Kmax *1.0/Settings.AccumulationFactor,2);
	// }
	// return base;
	double base = -1e3*pow(Parameters.Nu-83.0/Settings.Ploidy,2);
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

void Model::ComputeGradient(StateVector & grad,const std::vector<int> & histogram)
{
	grad.Reset();
	double epsInv = 1.0 - Parameters.Epsilon;
	double epsFactor = Parameters.Epsilon/(1.0 - Parameters.Epsilon) * (1.0 - Parameters.Epsilon/Settings.ErrorMax);
	int prevRhoIndex = -1;
	double expERho = 0;
	std::fill(digammaArray.begin(),digammaArray.end(),0.0);
	double dmin = 0;
	double dmax = Settings.ContaminationMax;
	for (int q = 0; q < NHarmonic; ++q)
	{
		contaminationPref[q] = (Parameters.Contamination[q] - dmin) * (dmax - Parameters.Contamination[q])/((dmax - dmin) * (q + Parameters.Contamination[q])); 

		dmin = Settings.ContaminationMin;
	}
	for (int k = 0; k <= Kmax; ++k)
	{
		if (histogram[k] > 0)
		{
			
			
			double dLdp = exp(log(histogram[k]) - ProbabilityArray[k]);
			double wSum = 0;
			for (int q = 0; q < NHarmonic; ++q)
			{
				wSum += Parameters.Weight[q] * exp(logB[q][k]);
			}
			
			
			for (int q = 0; q < NHarmonic; ++q)
			{
				grad.z[q] += dLdp * (1.0 - Parameters.Epsilon) * Parameters.Weight[q] * (exp(logB[q][k]) - wSum);

		
				double tau = Rq[q] * (digammaArray[q] - digammaMean[q]);
				double lambda = Pq[q]* (sMean[q]- k);
				
				
				double bcont = epsInv * exp(Parameters.LogWeight[q] + logB[q][k]);
				grad.x += bcont * (2 * tau + lambda) * dLdp;
				grad.y -= bcont * (tau + lambda) * dLdp;
				// LOG(DEBUG) << std::setprecision(8) <<"  " << k << " " << q << " " << tau << " " << lambda << "  " << dLdp << " " << sMean[q] << " " << grad.x;
				if (q != Settings.Ploidy)
				{
					grad.psi[q] += bcont * (2 * tau + lambda) * contaminationPref[q] * dLdp;
				}
			}
			// exit(5);
			//Error Quanitty
			grad.phi += epsFactor * (exp(LogError(k)) - exp(ProbabilityArray[k])) * dLdp;

			//Error Shape
			int RhoIndex = ErrorWindow(k);
			if (RhoIndex != prevRhoIndex)
			{
				expERho = exp(Parameters.E[RhoIndex]);
			}
			for (int j = 0; j < Settings.ErrorRes; ++j)
			{
				double mod = -exp(Parameters.E[j]);
				double width = Parameters.StandardWidth;
				if (j == Settings.ErrorRes-1)
				{
					width = Parameters.ResidualWidth;
				}
				mod *= width;
				if (j == RhoIndex)
				{
					mod +=1;
				}

				grad.gamma[j] += Parameters.Epsilon * expERho * mod * dLdp;
			}

		}

		for (int q = 0; q < NHarmonic; ++q)
		{
			digammaArray[q] += 1.0/(k + Rq[q]);
		}
	}
	
	// LOG(WARN) << JSL::Vector(grad.gamma);
}