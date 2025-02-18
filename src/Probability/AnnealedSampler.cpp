#include "AnnealedSampler.h"

AnnealedSampler::AnnealedSampler(const DataHolder & data)
{
	Histogram = data.Histogram();
	Vector = OptimiserParameters(7,Settings.ErrorRes);
	Proposed = OptimiserParameters(7,Settings.ErrorRes);
}

Model AnnealedSampler::Fit()
{
	Model P(Histogram.size()-1,Vector.z.size(),Settings.AccumulationFactor,Settings.ErrorRes);

	Random R;
	OptimiserParameters best;
	int Nit = 2000;
	// std::vector<double> mu;
	// std::vector<double> scores;

	P.SetParameters(Vector);
	double prevScore = abs(P.Score(Histogram));
	double minScore = prevScore;
	int rejectRun = 0;
	double T = Nit + 10;
	int overheats = 0;
	double stepSize = 0.5;
	double coolRate = pow(T,-1.0/Nit);
	int nansEncountered = 0;
	for (int i = 0; i < Nit; ++i)
	{
		// Vector.x = xs[i];
		Vector.RandomStep(R,Proposed,stepSize);
		P.SetParameters(Proposed);
		// mu.push_back(P.Parameters.Nu);
		double s = abs(P.Score(Histogram));
		// scores.push_back(s);
		double r = R.UniformDouble(0,1);

		
		double E = (prevScore - s) / T;
		if (std::isnan(E))
		{
			++nansEncountered;
			if (nansEncountered < 5)
			{
				LOG(WARN) << "NAN Encountered in optimisation routine. Attempting to recover.";
				Vector = best;
				prevScore = minScore;
			}
			else
			{
				throw runtime_error("Encountered 5 or more NAN errors. This signifies a major error.");
			}
		}
		else
		{
			if (r < exp(min(0.0,E)))
			{
				Vector = Proposed;
				prevScore = s;
				rejectRun = 0;
				overheats = 0;
			}
			else
			{
				++rejectRun;
				if (rejectRun >= 40)
				{
					// Vector = best;
					LOG(INFO) << i <<" Reheating from T=" << T << " " << stepSize;
					T *= 100;
					rejectRun =0;
					++overheats;
					stepSize = 0.1;
					coolRate =  pow(T,-1.0/(Nit-i));
					// stepSize = min(0.6,stepSize);
				}
				if (overheats > 1)
				{
					LOG(WARN) << "Overheat risk!";
					T *= 10;
					Vector = best;
					prevScore = minScore;
					rejectRun = 0;
					overheats = 0;
					// stepSize = 0.5;
				}
			}
			
			if (s < minScore)
			{
				minScore = s;
				best = Vector;
				LOG(DEBUG) << "New best " << minScore << " " << P.Parameters.Nu;
			}
		}
		T *= coolRate;
		stepSize *= 0.999;
	}

	// JSL::gnuplot gp;
	// gp.Scatter(mu,scores);
	// gp.SetYLog(true);
	// gp.SetXLog(true);
	// gp.Show();

	Vector = best;
	P.SetParameters(Vector);
	LOG(WARN) << "best x at " << P.Parameters.Nu;
	std::vector<double> ks = JSL::Vector::intspace(0,Histogram.size()-1,1);
	std::vector<int> plotx = JSL::Vector::intspace(0,ks.size()/Settings.AccumulationFactor,1);
	std::vector<double> mod(plotx.size(),0.0);
	std::vector<double> err(plotx.size(),0.0);
	std::vector<double> plotFreq(plotx.size(),0.0);
	int N = JSL::Vector(Histogram).Sum();
	double p = 0;
	double s = 0;
	for (int k: ks)
	{
		int roundk = k/Settings.AccumulationFactor;
		mod[roundk] +=exp(P[k] + log(N));
		plotFreq[roundk] += Histogram[k];
		err[roundk] += exp(P.LogError(k) + log(N)) * P.Parameters.Epsilon;
		// s = ale(s,P[k] -P.LogError(k));
		p += exp(P[k]);
		s += exp(P.LogError(k));
	}
	LOG(ERROR) << p << " " << s << " " << s/(p);
	LOG(WARN) << P.Parameters.Nu << " " << P.Parameters.Variance << " " << P.Parameters.Epsilon;
	LOG(WARN) << JSL::Vector(P.Parameters.Weight);
	LOG(WARN) << P.Parameters.Eta << "  " << JSL::Vector(P.Parameters.Contamination);
	JSL::gnuplot gp2;
	// // LOG(DEBUG) << ks.size() << " " << Histogram.size();
	gp2.Plot(plotx,plotFreq);
	gp2.Plot(plotx,mod);
	gp2.Plot(plotx,err,JSL::LineProperties::PenType(JSL::Dash));
	gp2.SetYLog(true);
	// gp2.SetXLog(true);
	gp2.Show();

	return P;
}