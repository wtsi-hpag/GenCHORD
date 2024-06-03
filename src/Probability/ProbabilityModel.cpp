#include "ProbabilityModel.h"


ProbabilityModel::ProbabilityModel(int dimensions)
{
	Dimensionality = dimensions;
}

double ProbabilityModel::logNoise(int k)
{
	return -99999;
}

double ProbabilityModel::logSignal(int k, int q)
{
	return -99999;
}


void ProbabilityModel::SetNoiseParameters_Observed(double weight, double mean, double sigma)
{
	NoiseWeight = weight;

	auto ms = GetMuSigma(mean,sigma);
	NoiseMean = ms.Mu; //GetMu(mean,sigma);
	NoiseSigma = ms.Sigma;//GetSigma(mini,minj,sigma);
}

MuSigmaPair ProbabilityModel::GetMuSigma(double muObs, double sigmaObs)
{
	double dmu = (inferredMeanMax-inferredMeanMin)/NoiseResolution;
	double dsig = (inferredSigmaMax-inferredSigmaMin)/NoiseResolution;
	
	muObs = std::min(inferredMeanMax-dmu,std::max(inferredMeanMin,muObs));
	sigmaObs = std::min(inferredSigmaMax-dsig,std::max(inferredSigmaMin,sigmaObs));

	int mu_coord = (muObs - inferredMeanMin)/dmu;
	int sigma_coord = (sigmaObs - inferredSigmaMin)/dsig;

	if (mu_coord == NoiseResolution -1)
	{
		--mu_coord;
	}
	if (sigma_coord == NoiseResolution-1)
	{
		--sigma_coord;
	}
	double muInterp = (muObs - inferredMeanMin - mu_coord*dmu)/dmu;	
	double sigmaInterp = (sigmaObs -inferredSigmaMin - sigma_coord*dsig)/dsig;

	auto m1 =  NoiseGrid[mu_coord][sigma_coord];
	auto m2 =  NoiseGrid[mu_coord+1][sigma_coord];
	auto m3 = NoiseGrid[mu_coord][sigma_coord+1];
	auto m4 = NoiseGrid[mu_coord+1][sigma_coord+1];

	double muLower = m1.Mu + (m2.Mu - m1.Mu) *muInterp;
	double sigmaLower = m1.Sigma + (m2.Sigma - m2.Sigma)*muInterp;

	double muUpper = m3.Mu + (m4.Mu - m3.Mu)*muInterp;
	double sigmaUpper = m3.Sigma + (m4.Sigma - m4.Sigma)*muInterp;

	double mu = muLower + (muUpper - muLower) * sigmaInterp;
	double sigma = sigmaLower + (sigmaUpper - sigmaLower)*sigmaInterp;

	// std::cout << "\t\tI think mean should be " << muObs << " I think this is given by " << m1.Mu << " " << m2.Mu << " " << m3.Mu << " " << m4.Mu << " " << mu << std::endl;
	// std::cout << "\t\tI think var should be " << sigmaObs << " I think this is given by " << m1.Sigma << " " << m2.Sigma << " " << m3.Sigma << " " << m4.Sigma << " " << sigma << std::endl;
	return MuSigmaPair(mu,sigma);

}



MuSigmaPair ProbabilityModel::ComputeMoments(int kMax)
{
	double s = 0;
	double sq = 0;

	for (int k = 0; k <= kMax; ++k)
	{
		double p = exp(logNoise(k) - NoiseNormalisation);
		s += k*p;
		sq += k*k*p;
	}
	
	return MuSigmaPair(s,sqrt(sq-s*s));
}

void ProbabilityModel::ComputeNoiseConversionChart(int kMax)
{
	double origMu = NoiseMean;
	double origSigma = NoiseSigma;
	std::vector<double> spoofws(0);


	std::vector<std::vector<MuSigmaPair>> SampleGrid(NoiseResolution);
	NoiseGrid.resize(NoiseResolution);


	double minMu = 100000;
	double maxMu = -1;
	double minSig = 1e10;
	double maxSig = -1;
	for (int nm = 0; nm < NoiseResolution; ++nm) 
	{
		SampleGrid[nm].resize(NoiseResolution);
		NoiseGrid[nm].resize(NoiseResolution);
		double trueMu = NoiseMeanMin + nm*(NoiseMeanMax - NoiseMeanMin)/(NoiseResolution - 1);

		for (int ns = 0; ns < NoiseResolution; ++ns)
		{
			double trueSigma = NoiseSigmaMin + ns * (NoiseSigmaMax - NoiseSigmaMin)/(NoiseResolution -1);

			//iniitial guesses
			NoiseMean = trueMu;
			NoiseSigma = trueSigma;
			Normalise(kMax,spoofws);
			
			SampleGrid[nm][ns] = ComputeMoments(kMax);

			double mu = SampleGrid[nm][ns].Mu;
			double sigma = SampleGrid[nm][ns].Sigma;

			if (mu > maxMu)
			{
				maxMu = mu;
			}
			if (mu < minMu)
			{
				minMu = mu;
			}
			if (sigma > maxSig)
			{
				maxSig = sigma;
			}
			if (sigma < minSig)
			{
				minSig = sigma;
			}
		}	
	}

	inferredMeanMin = minMu;
	inferredMeanMax = maxMu;
	inferredSigmaMin = minSig;
	inferredSigmaMax = maxSig;

	JSL::ProgressBar<2> PB(NoiseResolution,NoiseResolution);
	Log("Initialising truncated-distribution lookup table\n");
	for (int nm = 0; nm < NoiseResolution; ++nm) 
	{
		double observedMu = minMu + nm*(maxMu - minMu)/(NoiseResolution - 1);

		for (int ns = 0; ns < NoiseResolution; ++ns)
		{
			double observedSigma = minSig + ns*(maxSig - minSig)/(NoiseResolution - 1);

			double mind = 1e10;
			int mini;
			int minj;
			for (int i = 0; i < NoiseResolution; ++i)
			{
				for (int j = 0; j < NoiseResolution; ++j)
				{
					double d1 = (SampleGrid[i][j].Mu - observedMu);
					double d2 = (SampleGrid[i][j].Sigma - observedSigma);

					double d = d1*d1 + d2*d2;
					if (d < mind)
					{
						mind = d;
						mini = i;
						minj = j;
					}
				}
			}


			// NoiseGrid[mini]

			if (mini == NoiseResolution-1)
			{
				--mini;
			}
			if (minj == NoiseResolution-1)
			{
				--minj;
			}


			double dmu = (NoiseMeanMax-NoiseMeanMin)/(NoiseResolution-1);
			double dsig = (NoiseSigmaMax- NoiseSigmaMin)/(NoiseResolution-1);

			double muZero = NoiseMeanMin + mini*dmu;
			double sigZero = NoiseSigmaMin + minj*dsig;

			double muInterp = (SampleGrid[mini][minj].Mu - observedMu)/ (SampleGrid[mini+1][minj].Mu - (SampleGrid[mini][minj].Mu - muZero));
			double sigmaInterp = (SampleGrid[mini][minj].Sigma - observedSigma)/(SampleGrid[mini][minj+1].Sigma - (SampleGrid[mini][minj].Sigma - muZero));
			// double 
			muInterp = std::min(1.2,std::max(muInterp,-0.2));
			sigmaInterp= std::min(1.2,std::max(sigmaInterp,-0.2));

			double predictMu = std::min(NoiseMeanMax,std::max(NoiseMeanMin,muZero + dmu * muInterp));
			double predictSigma = std::min(NoiseSigmaMax,std::max(NoiseSigmaMin,sigZero + dsig * sigmaInterp));
			NoiseGrid[nm][ns] = MuSigmaPair(predictMu,predictSigma);

			PB.Update(nm,ns);
		}
	}

	NoiseMean = origMu;
	NoiseSigma = origSigma;
}



double ProbabilityModel::logP(int k, int q)
{
	return logP(k,q,true);
}
double ProbabilityModel::logP(int k, int q,bool withCorrections)
{
	double noiseWeighting = log(std::max(NoiseWeight,1e-100));
	double signalWeighting = log(1.0 - NoiseWeight);
	double signal = logSignal(k,q);
	double noise = logNoise(k);

	if (withCorrections)
	{
		signal -= Normalisation[q];
		noise -= NoiseNormalisation;
	}
	return ale(noiseWeighting+noise, signalWeighting+signal);
}

void ProbabilityModel::SetParameters(std::vector<double> v)
{
	Parameters = v;
}


void ProbabilityModel::Normalise(int kMax, std::vector<double> & weights)
{
	if (Normalisation.size() != weights.size())
	{
		Normalisation = std::vector<double>(weights.size(),0.0);
	}

	NoiseNormalisation = -999999;
	for (int q = 0; q < Normalisation.size(); ++q)
	{
		Normalisation[q] = -9999999;
	}
	for (int k = 0; k <= kMax; ++k)
	{
		NoiseNormalisation = ale(NoiseNormalisation,logNoise(k));
		for (int q = 0; q < Normalisation.size(); ++q)
		{
			Normalisation[q] = ale(Normalisation[q],logSignal(k,q));
		}
	}

}

void ProbabilityModel::FillGrid(std::vector<std::vector<double>> & grid)
{
	double noiseWeighting = log(std::max(NoiseWeight,1e-100));
	double signalWeighting = log(1.0 - NoiseWeight);
	for (int k = 0; k < grid.size(); ++k)
	{
		for (int q = 0; q < grid[k].size(); ++q)
		{
			grid[k][q] = logP(k,q);
		}
	}
}

std::vector<double> ProbabilityModel::GlobalPrediction(int kMax, std::vector<double> & weights)
{
	std::vector<double> out(kMax+1,0.0);
	double noiseWeighting = log(std::max(NoiseWeight,1e-100));
	double signalWeighting = log(1.0 - NoiseWeight);
	for (int q = 0; q < weights.size(); ++q)
	{
		for (int k = 0; k <= kMax; ++k)
		{
			double logp = logP(k,q);
			out[k] += weights[q] * exp(logp);
		}
	}
	return out;
}