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
	int Nit = 12000;
	std::vector<double> mu;
	std::vector<double> scores;
	std::vector<double> acceptance;
	std::vector<double> Ts;
	std::vector<double> baseline;
	P.SetParameters(Vector);
	double prevScore = abs(P.Score(Histogram));
	double minScore = prevScore;
	int rejectRun = 0;
	double T = Nit + 10;
	int overheats = 0;
	double stepSize = 0.5;
	double finalStepSize = 0.05;
	double coolRate = pow(T,-2.0/Nit);
	double stepRate = pow(finalStepSize/stepSize,1.0/Nit);
	int nansEncountered = 0;
	double accept = 0;
	double reject = 0;
	int timeSinceHeat = 0;
	int timeSinceQuench = 0;
	double mem = (1.0 - 100/Nit);
	for (int i = 0; i < Nit; ++i)
	{
		baseline.push_back(i);
		// Vector.x = xs[i];
		Vector.RandomStep(R,Proposed,stepSize);
		P.SetParameters(Proposed);
		mu.push_back(P.Parameters.Nu);
		double s = abs(P.Score(Histogram));
		scores.push_back(s);
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
				accept = mem*(accept + 1);
				// coolRate =  pow(T,-1.0/(Nit-i));
				timeSinceHeat = max(0,timeSinceHeat -1);
				timeSinceQuench = max(0,timeSinceQuench-1);
			}
			else
			{
				reject = mem*(reject + 1);
				++rejectRun;
				if (accept/(accept + reject ) < 0.25 && timeSinceHeat > 20)
				{
					// Vector = best;
					T *= 10;
					rejectRun =0;
					++overheats;
					stepSize = 0.1;
					stepRate = pow(finalStepSize/stepSize,2.0/(Nit-i));
					coolRate =  pow(T,-1.0/(Nit-i));
					timeSinceHeat =0;
					// stepSize = min(0.6,stepSize);
				}
				if (overheats > 1)
				{
					T *= 10;
					Vector = best;
					prevScore = minScore;
					rejectRun = 0;
					overheats = 0;
				}
			}
			
			if (s < minScore)
			{
				minScore = s;
				best = Vector;
			}
		}
		++timeSinceHeat;
		++timeSinceQuench;
		Ts.push_back(T);
		double moverate = accept/(accept + reject);
		acceptance.push_back(moverate);
		if (moverate > 0.5 && timeSinceQuench > 30)
		{
			coolRate = max(0.8,coolRate * 0.99);
			timeSinceQuench = 0;
		}
		
		T = max(1.0,T*coolRate);
		if (T == 1)
		{
			coolRate = 1.01;
		}
		stepSize *= stepRate;
	}
	LOG(WARN) << "Accepted = " << 100 * accept *1.0/ (accept + reject);

	for (auto & x: scores)
	{
		x = 1 + x - minScore;
	}

	JSL::gnuplot gp;
	gp.SetMultiplot(3,1);
	gp.WindowSize(800,600);
	gp.Scatter(mu,scores);
	gp.SetYLog(true);
	gp.SetXLog(true);
	gp.SetAxis(1);
	gp.Plot(baseline,acceptance);
	gp.SetTitle("Acceptance Fraction");
	gp.SetAxis(2);
	gp.Plot(baseline,Ts);
	gp.SetTitle("Temperature");
	gp.SetYLog(true);
	gp.Show();

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
	int maxObs = 0;
	for (int k: ks)
	{
		int roundk = k/Settings.AccumulationFactor;
		mod[roundk] +=exp(P[k] + log(N));
		plotFreq[roundk] += Histogram[k];
		err[roundk] += exp(P.LogError(k) + log(N)) * P.Parameters.Epsilon;
		// s = ale(s,P[k] -P.LogError(k));
		p += exp(P[k]);
		s += exp(P.LogError(k));
		if (plotFreq[roundk] > maxObs)
		{
			maxObs = plotFreq[roundk];
		}
		if (mod[roundk] > maxObs)
		{
			maxObs = mod[roundk];
		}
	}
	LOG(ERROR) << p << " " << s << " " << s;
	LOG(WARN) << P.Parameters.Nu << " " << P.Parameters.Variance << " " << P.Parameters.Epsilon;
	LOG(WARN) << JSL::Vector(P.Parameters.Weight);
	LOG(WARN) << P.Parameters.Eta << "  " << JSL::Vector(P.Parameters.Contamination);
	JSL::gnuplot gp2;
	// // LOG(DEBUG) << ks.size() << " " << Histogram.size();
	gp2.Plot(plotx,plotFreq);
	gp2.Plot(plotx,mod);
	gp2.Plot(plotx,err,JSL::LineProperties::PenType(JSL::Dash));

	for (int q = 0; q < P.NHarmonic; ++q)
	{
		double p =0;
		std::vector<double> harmonic(plotx.size());
		for (int k: ks)
		{
			int roundk = k/Settings.AccumulationFactor;
			harmonic[roundk] +=exp(P.Sample(q,k) + log(N));
			p += exp(P.Sample(q,k));
		}
		LOG(WARN) << "q = " << q << " normed to " << p/((1 - P.Parameters.Epsilon) * P.Parameters.Weight[q]);
		gp2.Plot(plotx,harmonic,JSL::LineProperties::Legend("q = " + std::to_string(q)),JSL::LineProperties::PenType(JSL::Dotted));
	}
	gp2.SetYRange(0.1,maxObs);
	gp2.SetLegend(true);
	gp2.SetYLog(true);
	// gp2.SetXLog(true);
	gp2.Show();

	return P;
}