#include "ProbabilityModel.h"

ProbabilityModel::ProbabilityModel()
{

}
ProbabilityModel::ProbabilityModel(int kmax, int qmax)
{
	SetDimensionality(kmax,qmax);
}
void ProbabilityModel::SetDimensionality(int kmax, int qmax)
{
	MaxK = kmax;
	MaxQ = qmax;

	Noise = std::vector<double>(kmax+1,0.0);
	Signal = std::vector<std::vector<double>>(kmax+1,std::vector<double>(qmax+1,0.0));
}
Dual ProbabilityModel::GetDimensionality()
{
	return Dual(MaxK,MaxQ);
}

double ProbabilityModel::logNoise(int k)
{
	return -99999;
}

double ProbabilityModel::logSignal(int k, int q)
{
	return -99999;
}

void ProbabilityModel::SetSignalParameters(double mu, double sigma)
{
	SignalMean = mu;
	SignalSigma = sigma;
}
void ProbabilityModel::SetSignalParameters(Dual input)
{
	SetSignalParameters(input.X,input.Y);
}
void ProbabilityModel::SetNoiseParameters(double mu, double sigma, double gamma)
{

	NoiseMean = mu;
	NoiseSigma = sigma;
	NoiseWeight = std::max(1e-100,gamma);
	logNoiseWeight = log(NoiseWeight);
	logSignalWeight = log(1.0 - NoiseWeight);
}

void ProbabilityModel::SetGrids()
{
	SetNoiseGrid();
	SetSignalGrids();
}
void ProbabilityModel::SetNoiseGrid()
{
	double noiseNorm = -999999999;
	for (int k = 0; k <=MaxK; ++k)
	{
		Noise[k] = logNoise(k);
		noiseNorm = ale(noiseNorm,Noise[k]);
	}
	for (int k = 0; k <=MaxK; ++k)
	{
		Noise[k] -= noiseNorm;
	}
}
void ProbabilityModel::SetSignalGrids()
{
	for (int q = 0; q<=MaxQ; ++q)
	{
		double noiseNorm = -99999999999;
		for (int k = 0; k <=MaxK; ++k)
		{
			Signal[k][q] = logSignal(k,q);
			noiseNorm = ale(noiseNorm,Signal[k][q]);
		}
		for (int k = 0; k <=MaxK; ++k)
		{
			Signal[k][q] -= noiseNorm;
		}
	}
}



void ProbabilityModel::SetNoiseParametersFromObserved(double mean, double sigma, double weight)
{
	auto ms = GetMuSigma(mean,sigma);
	SetNoiseParameters(ms.X,ms.Y,weight);
}

Dual ProbabilityModel::GetMuSigma(double muObs, double sigmaObs)
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

	double muLower = m1.X + (m2.X - m1.X) *muInterp;
	double sigmaLower = m1.Y + (m2.Y - m2.Y)*muInterp;

	double muUpper = m3.X + (m4.X - m3.X)*muInterp;
	double sigmaUpper = m3.Y + (m4.Y - m4.Y)*muInterp;

	double mu = muLower + (muUpper - muLower) * sigmaInterp;
	double sigma = sigmaLower + (sigmaUpper - sigmaLower)*sigmaInterp;
	return Dual(mu,sigma);

}



Dual ProbabilityModel::ComputeNoiseMoments()
{
	double s = 0;
	double sq = 0;

	for (int k = 0; k <= MaxK; ++k)
	{
		double p = exp(Noise[k]);
		s += k*p;
		sq += k*k*p;
	}
	
	return Dual(s,sqrt(sq-s*s));
}

void ProbabilityModel::ComputeNoiseConversionChart()
{
	double origMu = NoiseMean;
	double origSigma = NoiseSigma;


	std::vector<std::vector<Dual>> SampleGrid(NoiseResolution);
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
			SetNoiseGrid();
			
			SampleGrid[nm][ns] = ComputeNoiseMoments();

			double mu = SampleGrid[nm][ns].X;
			double sigma = SampleGrid[nm][ns].Y;

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
					double d1 = (SampleGrid[i][j].X - observedMu);
					double d2 = (SampleGrid[i][j].Y - observedSigma);

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

			double muInterp = (SampleGrid[mini][minj].X - observedMu)/ (SampleGrid[mini+1][minj].X - (SampleGrid[mini][minj].X - muZero));
			double sigmaInterp = (SampleGrid[mini][minj].Y - observedSigma)/(SampleGrid[mini][minj+1].Y - (SampleGrid[mini][minj].Y - muZero));
			// double 
			muInterp = std::min(1.2,std::max(muInterp,-0.2));
			sigmaInterp= std::min(1.2,std::max(sigmaInterp,-0.2));

			double predictMu = std::min(NoiseMeanMax,std::max(NoiseMeanMin,muZero + dmu * muInterp));
			double predictSigma = std::min(NoiseSigmaMax,std::max(NoiseSigmaMin,sigZero + dsig * sigmaInterp));
			NoiseGrid[nm][ns] = Dual(predictMu,predictSigma);

			PB.Update(nm,ns);
		}
	}

	NoiseMean = origMu;
	NoiseSigma = origSigma;
}


double ProbabilityModel::getLogNoise(int k)
{
	return Noise[k];
}
double ProbabilityModel::getLogSignal(int k, int q)
{
	return Signal[k][q];
}
double ProbabilityModel::logP(int k, int q)
{
	// return Noise[k];
	return ale(logNoiseWeight + Noise[k], logSignalWeight + Signal[k][q]);
}

void ProbabilityModel::GlobalLogPrediction(std::vector<double> & out, const std::vector<double> & weights)
{
	if (out.size() != MaxK+1)
	{
		std::cout << "ERROR! Kmax of probability and output do not match. Ensuring Kmax is set properly is crucial for probability normalisation." << std::endl;
		exit(5);
	}

	for (int k = 0; k <= MaxK; ++k)
	{
		for (int q = 0; q < weights.size(); ++q)
		{
		
			double p_kq = log(weights[q]) +  logP(k,q);
			if (q == 0)
			{
				out[k] = p_kq;
			}
			else
			{
				out[k] = ale(p_kq,out[k]);
			}
		}
	}
}