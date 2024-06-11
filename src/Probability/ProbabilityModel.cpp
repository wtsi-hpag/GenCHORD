#include "ProbabilityModel.h"
#include "../../libs/JSL/JSL.h"
ProbabilityModel::ProbabilityModel(int noiseResolution)
{
	NoiseResolution = noiseResolution;
}
ProbabilityModel::ProbabilityModel(int noiseResolution,int kmax, int qmax)
{
	NoiseResolution = noiseResolution;
	SetDimensionality(kmax,qmax);
}
void ProbabilityModel::SetDimensionality(int kmax, int qmax)
{
	MaxK = kmax;
	MaxQ = qmax;

	if (Contamination.size() != qmax+1)
	{
		Contamination = std::vector<double>(qmax+1,0.0);
	}
	Noise = std::vector<double>(kmax+1,0.0);
	Signal = std::vector<std::vector<double>>(kmax+1,std::vector<double>(qmax+1,0.0));
	NoiseComponents = std::vector<double>(NoiseResolution);
	logNoiseWidth = std::vector<double>(NoiseResolution,0.);
	NoiseResolution = std::min(MaxK,NoiseResolution);
	int prev = 0;
	for (int r = 0; r < NoiseResolution; ++r)
	{
		int typical = MaxK/NoiseResolution;
		int insert = typical;
		if (r == NoiseResolution - 1)
		{
			insert = MaxK+1 - (NoiseResolution-1)*typical;
		}
		logNoiseWidth[r] = log(insert);
		prev += insert;
	}
	ResetNoise();
}
Dual ProbabilityModel::GetDimensionality()
{
	return Dual(MaxK,MaxQ);
}
void ProbabilityModel::ResetNoise()
{

	double rSum =0;
	for (int r = 0; r < NoiseResolution; ++r)
	{
		NoiseComponents[r] = 1;
		rSum += NoiseComponents[r];
	}
	for (int r = 0; r < NoiseResolution; ++r)
	{
		NoiseComponents[r] = log(NoiseComponents[r]/rSum);
	}
	logNoiseWeight = -4;
	logSignalWeight = log(1.0 - exp(logNoiseWeight));
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

int ProbabilityModel::RFromK(int k)
{
	int ksPerBlock = MaxK * 1.0/NoiseResolution;
	return std::min(NoiseResolution-1,k/ksPerBlock);
}
void ProbabilityModel::SetGrids()
{
	SetNoiseGrid();
	SetSignalGrids();
	// SetConvolutionGrids();
}
void ProbabilityModel::SetNoiseGrid()
{
	// double noiseNorm = -999999999;
	for (int k = 0; k <=MaxK; ++k)
	{
		Noise[k] = logNoise(k);
		// noiseNorm = ale(noiseNorm,Noise[k]);
	}
	// for (int k = 0; k <=MaxK; ++k)
	// {
	// 	Noise[k] -= noiseNorm;
	// }
}
void ProbabilityModel::SetSignalGrids()
{
	for (int q = 0; q<=MaxQ; ++q)
	{
		double noiseNorm = -99999999999;
		for (int k = 0; k <=MaxK; ++k)
		{
			int kCorrect = k;
			// if (q != 2 && q != 0)
			// {
			// 	kCorrect = std::max(0,(int)round(k -Contamination[q]*2*SignalMean));
			// }
			Signal[k][q] = logSignal(kCorrect,q);
			noiseNorm = ale(noiseNorm,Signal[k][q]);
		}
		for (int k = 0; k <=MaxK; ++k)
		{
			Signal[k][q] -= noiseNorm;
		}
	}
}
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
void ProbabilityModel::SetContamination(int q, double targetMean)
{
	double old = targetMean;
	double maxShift = 0.5*SignalMean;
	if (q == 0)
	{
		maxShift = 5;
	}
	double shift = std::min(maxShift,abs(q*SignalMean - targetMean));
	double amount = shift *sgn(-q*SignalMean +targetMean);
	targetMean = amount + q*SignalMean;
	
	// std::cout << q*SignalMean << "  " << old << "  " << shift << "  " << amount << "  " << targetMean << std::endl;
	double impliedCont = (targetMean - q*SignalMean)/((2-q)*SignalMean);
	// impliedCont = std::max(0.,std::min(0.5,impliedCont));
	Contamination[q] = impliedCont;
	// if (q == 0)
	// {
	// 	Contamination[q] = targetMean;
	// }
	// else
	// {
	// 	if (targetMean > currentMean)
	// 	{
	// 		double impliedCont = (targetMean - currentMean)/(2*SignalMean);
	// 		if (impliedCont < 0.1)
	// 		{
	// 			Contamination[q] = 0.5*Contamination[q] + 0.5*impliedCont;
	// 		}
	// 	}
	// }
}



void ProbabilityModel::ClearContamination()
{
	std::fill(Contamination.begin(),Contamination.end(),0.);
	// Contamination[1] = 0.05;
	// Contamination[3] = 0.4;
	// Contamination[0] = 0.01;
}


double ProbabilityModel::getLogNoise(int k)
{
	if (k > MaxK)
	{
		k = MaxK;
	}
	// std::cout << "Noise " << k << "  " << Noise[k] << std::endl;
	return Noise[k];
}
double ProbabilityModel::getLogNoise(int k, int r)
{
	if (k > MaxK)
	{
		k = MaxK;
	}
	int ksPerBlock = MaxK * 1.0/NoiseResolution;
	int kBelongsTo = std::min(NoiseResolution-1,k/ksPerBlock);
	if (kBelongsTo!=r)
	{
		return -9999999999;
	}
	else
	{
		return Noise[k];
	}
}


double ProbabilityModel::getLogSignal(int k, int q)
{
	if (k > MaxK)
	{
		k = MaxK;
	}
	return Signal[k][q];
}
double ProbabilityModel::logP(int k, int q)
{
	// return Noise[k];
	if (k > MaxK)
	{
		k = MaxK;
	}
	return ale(logNoiseWeight+Noise[k], logSignalWeight + Signal[k][q]);
}

void ProbabilityModel::SetNoiseComponents(std::vector<double> & wr)
{
	double wSum = 0;
	for (int r = 0; r < wr.size(); ++r)
	{
		wSum += wr[r];
	}
	for (int r = 0; r < wr.size(); ++r)
	{
		NoiseComponents[r] = log(wr[r]/wSum);
	}
	// double gamma = 0;
	// for (int r =0 ; r < NoiseComponents.size(); ++r)
	// {
	// 	NoiseComponents[r] = log(std::max(1e-100,wr[r]));
	// 	gamma += wr[r];
	// }
	// double gammaCrit = 0.05;
	// if (gamma > gammaCrit)
	// {
	// 	for (int r = 0; r < NoiseComponents.size(); ++r)
	// 	{
	// 		NoiseComponents[r] -= log(gamma/gammaCrit);
	// 	}
	// 	gamma = gammaCrit;
	// }
	// std::cout << "new gamma " << gamma << std::endl;
	// logSignalWeight = log(1.0 - gamma);

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
			// std::cout << 
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