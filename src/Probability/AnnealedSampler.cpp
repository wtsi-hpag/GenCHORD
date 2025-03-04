#include "AnnealedSampler.h"


AnnealedSampler::AnnealedSampler(const DataHolder & data): Vector(Settings.HarmonicCount,Settings.ErrorRes)
{
	Histogram = data.Histogram();
	Proposed = StateVector(Settings.HarmonicCount,Settings.ErrorRes);
}

Model AnnealedSampler::Fit()
{
	Model P(Histogram.size()-1,Settings.HarmonicCount,Settings.AccumulationFactor,Settings.ErrorRes);

	Random R;
	StateVector best;
	int Nit = 12000;
	std::vector<double> mu;
	std::vector<double> scores;
	std::vector<double> acceptance;
	std::vector<double> Ts;
	std::vector<double> baseline;
	P.SetParameters(Vector);
	
	double prevScore = abs(P.Score(Histogram));
	double minScore = prevScore;
	double T = abs(prevScore);
	double stepSize = 0.1;
	double coolRate = pow(T,-2.0/Nit);
	int nansEncountered = 0;
	double accept = 0;
	int timeSinceHeat = 0;
	int timeSinceQuench = 0;
	double mem = 1.0 - 50.0/Nit;
	for (int i = 0; i < Nit; ++i)
	{
		baseline.push_back(i);
		Vector.Parameters.RandomStep(R,Proposed,stepSize);
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
				Vector.Parameters = best;
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
				//accept the proposal
				Vector.Parameters = Proposed;
				prevScore = s;
				accept = mem*accept + (1.0 - mem);

				//cool the system if too many positions accepted
				if (accept > 0.5 && timeSinceQuench > 30)
				{
					T/=2;
					timeSinceQuench = 0;
					// accept = 0.5;
				}
			}
			else
			{
				accept = mem * accept;

				if (timeSinceHeat > 20)
				{
					//if doing *really badly*, jump to the best score to get back in the mix.
					if (accept < 0.1)
					{
						Vector.Parameters = best;
						prevScore = minScore;
						T *= 10;
						// accept = 0.1;
					}
					else if (accept < 0.25)
					{
						double minFac = 2;
						double maxFac = 10;
						double heatFactor = maxFac - (maxFac - minFac)/Nit * i;
						T *= heatFactor;
						coolRate =  pow(T,-1.0/(Nit-i));
						timeSinceHeat =0;
						// accept = 0.25;
					}
				}
			}
			
			if (s < minScore)
			{
				minScore = s;
				best = Vector.Parameters;
			}
		}
		++timeSinceHeat;
		++timeSinceQuench;
		Ts.push_back(T);
		acceptance.push_back(accept);
		
		
		T = max(1.0,T*coolRate);
		if (T == 1)
		{
			coolRate = 1.01;
		}
	}
	LOG(WARN) << "Accepted = " << 100 * accept;

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