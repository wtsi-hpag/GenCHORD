#include "NegativeBinomial.h"

NegativeBinomial::NegativeBinomial(int kMax, int resolution,double maxMu)
{
	Kmax = kMax;
	PureData.resize(kMax+1,std::vector<double>(resolution,0.0));
	dMu = maxMu * 1.0/resolution;
	Resolution = resolution;
}
void NegativeBinomial::Populate(double sigma)
{
	// sigma = 7;
	if (sigma < 1e-7)
	{
		sigma = 1e-7;
	}
	
	Sigma = sigma;

	
	double var = sigma*sigma;
	for (int k = 0; k <= Kmax; ++k)
	{
		//do the zeroth component because it involves log(0) if done naively
		double zeroVal = -9999999999;
		if (k == 0)
		{
			zeroVal = 0;
		}
		PureData[k][0] = zeroVal;


		for (int i = 1; i < Resolution; ++i)
		{
			double mu = i * dMu;

			double n = mu*mu/(var);
			double p = mu/(mu + var);
			double firstTerm = lgamma(n + k) - lgamma(n) - LogFac(k);
			double secondTerm = n * log(p) + k * log(1.0 - p);

			PureData[k][i] = (firstTerm + secondTerm);
		}
	}
}


TruncatedGaussian::TruncatedGaussian(int kMax, int bounder)
{
	Truncator = bounder;
	KMax = kMax;
	PureData.resize(kMax+1,std::vector<double>(kMax + bounder+2,0.0));
}
void TruncatedGaussian::Populate(double gamma)
{
	for (int kobs = 0; kobs <= KMax; ++kobs)
	{
		double s = -999999999;
		for (int k = std::max(0,kobs-Truncator); k <= KMax + Truncator; ++k)
		{
			double d = (k - kobs)/gamma;
			double v = - 0.5 * d * d;
			PureData[kobs][k] = v;
			s = ale(s,v);
		}
		for (int k = std::max(0,kobs-Truncator); k <= KMax + Truncator; ++k)
		{
			PureData[kobs][k] -= s;
		}
	}
}

ErroredBinomial::ErroredBinomial(int kMax, int res, int bounder, int q,double muMax, int nWork) : NB(kMax+bounder,res,muMax), TG(kMax,bounder), CurrentBracket(q), Resolution(res), qMax(q)
{
	nWorkers = nWork;
	Data.resize(kMax+1,std::vector<double>(res,0.0));
	Kmax = kMax;
	Bounder = bounder;


}


// void ErroredBinomial::Worker(int id)
// {
// 	int delta = Resolution/(nWorkers);
// 	int resStart = id*delta;
// 	int resEnd = std::min(resStart + delta,Resolution);
// 	mtx.lock();
// 	masterStart = std::max(masterStart,resEnd);
// 	mtx.unlock();
// 	while (true)
// 	{
// 		mtx.lock();
// 		int q = WorkerStatus[id];
// 		mtx.unlock();
// 		if (q == -1)
// 		{
// 			break;
// 		}
// 		if (q == 1)
// 		{
// 			PopulateChunk(resStart, resEnd);
// 			mtx.lock();
// 			WorkerStatus[id] = 0;
// 			mtx.unlock();
// 		}
		
// 	}
// 	WorkerStatus[id] = 0;
// }

void ErroredBinomial::PopulateChunk(int resMin, int resMax)
{
	for (int i = resMin; i < resMax; ++i)
	{
		for (int kobs = 0; kobs <= Kmax; ++kobs)
		{
			double s = -9999999999;
			for (int k = std::max(0,kobs-Bounder); k <= Kmax + Bounder; ++k)
			{
				double v = TG.PureData[kobs][k] + NB.PureData[k][i];
				s = ale(s,v);
			}
			Data[kobs][i] = s;
		}
	}
}
void ErroredBinomial::Populate(double gamma,double sigma, Distributor & dist)
{
	Sigma = sigma;
	TG.Populate(gamma);
	NB.Populate(sigma);
	
	dist.UpdateEB(this);
	dist.Signal(1);

	int delta = Resolution/(dist.WorkerCount+1);
	int end = dist.WorkerCount * delta;
	masterStart = std::min(end,Resolution);
	PopulateChunk(masterStart,Resolution); //do chunk on main thread too
	dist.Gather();
}

void ErroredBinomial::SetBracket(double nu)
{

	for (int q = 0; q < qMax; ++q)
	{
		double mu = q*nu;
		
		double r = std::min((double)Resolution-1,mu/NB.dMu);
		CurrentBracket.Bottom[q] = floor(r);
		CurrentBracket.Top[q] = std::min(Resolution-1,(int)ceil(r));
		CurrentBracket.Interpolate[q] = r - CurrentBracket.Bottom[q];
	}
}

double ErroredBinomial::GetProb(int q, int k)
{
	double top = Data[k][CurrentBracket.Top[q]];
	double bottom = Data[k][CurrentBracket.Bottom[q]];
	return bottom + CurrentBracket.Interpolate[q] * (top-bottom);
	// return v;
}

// void ErroredBinomial::CloseThreads()
// {
// 	for (int i = 0; i < WorkerStatus.size(); ++i)
// 	{
// 		mtx.lock();
// 		WorkerStatus[i] = -1;
// 		mtx.unlock();
// 	}

// 	for (int i = 0; i < Threads.size(); ++i)
// 	{
// 		Threads[i].join();
// 	}
// }