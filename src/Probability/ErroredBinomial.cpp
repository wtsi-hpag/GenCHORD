#include "ErroredBinomial.h"

ErroredBinomial::ErroredBinomial(int kMax, int res, int bounder, int q,double muMax, int nWork) : NB(kMax+bounder,res,muMax), TG(kMax,bounder), CurrentBracket(q), Resolution(res), qMax(q)
{
	nWorkers = nWork;
	Data.resize(kMax+1,std::vector<double>(res,0.0));
	Kmax = kMax;
	Bounder = bounder;


}


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