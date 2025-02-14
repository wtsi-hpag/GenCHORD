#include "AnnealedSampler.h"

AnnealedSampler::AnnealedSampler(const DataHolder & data)
{
	Histogram = data.Histogram();

	Vector = OptimiserParameters(7);
}

Model AnnealedSampler::Fit()
{
	Model P(Histogram.size()-1,Vector.z.size(),Settings.AccumulationFactor);

	
	OptimiserParameters best;
	int Nit = 20000;
	std::vector<double> mu;
	std::vector<double> scores;
	double bm = 0;
	double min= 9e99;
	for (int i = 0; i < Nit; ++i)
	{
		// Vector.x = xs[i];
		Vector.Randomise();
		P.SetParameters(Vector);
		mu.push_back(P.Parameters.Nu);
		scores.push_back(abs(P.Score(Histogram)));
		if (scores[i] < min)
		{
			min = scores[i];
			bm = P.Parameters.Nu;
			best = Vector;
		}
	}

	JSL::gnuplot gp;
	gp.Scatter(mu,scores);
	// // gp.Plot(mus,scores);
	gp.SetYLog(true);
	gp.SetXLog(true);
	gp.Show();

	// Vector.x = bm;
	Vector = best;
	P.SetParameters(Vector);
	LOG(WARN) << "best x at" << bm;
	std::vector<double> ks = JSL::Vector::intspace(0,Histogram.size()-1,1);
	std::vector<int> plotx = JSL::Vector::intspace(0,ks.size()/Settings.AccumulationFactor,1);
	std::vector<double> mod(plotx.size(),0.0);
	std::vector<double> plotFreq(plotx.size(),0.0);
	int N = JSL::Vector(Histogram).Sum();
	double p = 0;
	for (int k: ks)
	{
		int roundk = k/Settings.AccumulationFactor;
		mod[roundk] +=exp(P[k] + log(N) - P.Normalisation);
		plotFreq[roundk] += Histogram[k];
		p += exp(P[k] - P.Normalisation);
		// ks[k] /= Settings.AccumulationFactor; 
		// Histogram[k] += 0.1;
	}

	LOG(WARN) << mod.size() << " " << plotx.size();
	LOG(WARN) << p << " <- probsum" << " " << P.Normalisation;
	LOG(WARN) << P.Parameters.Nu << " " << P.Parameters.Variance << " " << P.Parameters.Epsilon;
	LOG(WARN) << JSL::Vector(P.Parameters.Weight);
	LOG(WARN) << P.Parameters.Eta << "  " << JSL::Vector(P.Parameters.Contamination);
	JSL::gnuplot gp2;
	// LOG(DEBUG) << ks.size() << " " << Histogram.size();
	gp2.Plot(plotx,plotFreq);
	gp2.Plot(plotx,mod);
	gp2.SetYLog(true);
	// gp2.SetXLog(true);
	gp2.Show();

	return P;
}