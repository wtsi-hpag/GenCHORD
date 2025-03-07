#include "AnnealedSampler.h"


void PlotModel(JSL::gnuplot & gp, const std::vector<int> Histogram, Model & P, double & maxObs, bool plotData,std::string title)
{
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
		mod[roundk] += exp(P[k]);
		plotFreq[roundk] += Histogram[k]* 1.0/N;
		err[roundk] += exp(P.LogError(k)) * P.Parameters.Epsilon;
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
	
	if (plotData)
	{
		gp.Plot(plotx,plotFreq,JSL::LineProperties::Legend("Observed Data"));
	}
	gp.Plot(plotx,mod,JSL::LineProperties::Legend(title));
	// // gp2.Plot(plotx,err,JSL::LineProperties::PenType(JSL::Dash));

	// for (int q = 0; q < P.NHarmonic; ++q)
	// {
	// 	double p =0;
	// 	std::vector<double> harmonic(plotx.size());
	// 	for (int k: ks)
	// 	{
	// 		int roundk = k/Settings.AccumulationFactor;
	// 		harmonic[roundk] +=exp(P.Sample(q,k) + log(N));
	// 		p += exp(P.Sample(q,k));
	// 	}
	// 	gp2.Plot(plotx,harmonic,JSL::LineProperties::Legend("q = " + std::to_string(q)),JSL::LineProperties::PenType(JSL::Dotted));
	// }
	gp.SetYRange(0.1/N,1.0);
	gp.SetLegend(true);
	gp.SetYLog(true);

}

AnnealedSampler::AnnealedSampler(const DataHolder & data): Vector(Settings.HarmonicCount,Settings.ErrorRes), Proposed(Settings.HarmonicCount,Settings.ErrorRes), Data(data)
{
	Histogram = data.Histogram();
	// Model tempModel(Histogram.size()-1,Settings.HarmonicCount,Settings.AccumulationFactor,Settings.ErrorRes);
	// tempModel.Parameters.Nu = 45;
	// tempModel.Parameters.Variance = 10;
	// tempModel.Parameters.Weight[3] = 
	// auto q = StateVector(Settings.HarmonicCount,Settings.ErrorRes);
	// q.SetDefaultValues();
	
	// q.y = log(185);
	// q.x = log(19);
	// q.psi[0] = 0.3;
	// // q.psi[4] = -0.4;
	// // q.phi = -2;
	// tempModel.SetParameters(q);
	// // tempModel.Compute();
	// LOG(ERROR) << JSL::Vector(q.psi);
	// LOG(ERROR) << JSL::Vector(tempModel.Parameters.Contamination);
	// Random R;
	// for (int k = 0; k < Histogram.size(); ++k)
	// {
	// 	Histogram[k] = max(0,(int)(exp(14 + tempModel[k])));
	// }

	
}

Model AnnealedSampler::Fit()
{
	Model P(Histogram.size()-1,Settings.HarmonicCount,Settings.AccumulationFactor,Settings.ErrorRes);

	Random R;
	StateVector best(Settings.HarmonicCount,Settings.ErrorRes);
	best.SetDefaultValues();
	int Nit = 500;
	std::vector<double> mu;
	std::vector<double> scores;
	std::vector<double> acceptance;
	std::vector<double> Ts;
	std::vector<double> baseline;
	P.SetParameters(Vector);
	
	double prevScore = abs(P.Score(Histogram));
	double minScore = prevScore;
	double T = abs(prevScore);
	double stepSize = 0.05;
	double coolRate = pow(T,-2.0/Nit);
	int nansEncountered = 0;
	double accept = 0;
	int timeSinceHeat = 0;
	int timeSinceQuench = 0;
	int timeSinceOptim = 0;
	double mem = max(0.9,1.0 - 10.0/Nit);
	
	for (int i = 0; i < Nit; ++i)
	{
		baseline.push_back(i);

		Vector.Parameters.RandomStep(R,Proposed,stepSize);
		P.SetParameters(Proposed);
		if (timeSinceOptim < 100)
		{
			
			++timeSinceOptim;
		}
		else
		{
			timeSinceOptim = 0;
			Optimise(P,10 + (300*i)/Nit);
			
			// T*=100;
		}
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
				Vector.Parameters.LightCopy(Proposed);
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
					timeSinceHeat = 0;
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
						// accept = 0.25;
					}
				}
			}
			
			if (s < minScore)
			{
				minScore = s;
				best.LightCopy(Vector.Parameters);

				P.SetParameters(best);
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

	Vector.Parameters.LightCopy(best);
	P.SetParameters(best);
	P.Compute();
	// P.SetParameters(Vector);

	JSL::gnuplot gp2;
	double maxObs = 0;
	
	PlotModel(gp2,Histogram,P,maxObs,true,"Sampled");
	double s = P.Score(Histogram);
	Optimise(P,200);
	double s2 = P.Score(Histogram);
	PlotModel(gp2,Histogram,P,maxObs,false,"Optimised");
	LOG(DEBUG) << "Final improvement " << s << " -> " << s2;
	gp2.Show();
	return P;
}

void AnnealedSampler::Optimise(Model & model, int Nsteps)
{
	auto vCopy = Vector;
	double alpha = 0.1;
	double b1 = 0.5;
	double b2 = 0.9;
	Vector.OptimiseReset();
	for (int i = 0; i < Nsteps; ++i)
	{
		model.ComputeGradient(Vector.Gradient,Histogram);
		
		Vector.AccumulateGradient(b1,b2);
		Vector.ADAMUpdate(alpha,b1,b2,i+1);
		model.SetParameters(Vector);
		// q = (q + 1) % Settings.HarmonicCount;
		// LOG(DEBUG) << model.Parameters.Nu << "  " << JSL::Vector(Vector.Gradient.psi);
	}
	Proposed = Vector.Parameters;
	Vector = vCopy;
}

Model AnnealedSampler::FineTune(Model model, int chromosome)
{
	Histogram = Data.Histogram(chromosome);
	Data.TruncateHistogram(Histogram,Data[chromosome].size());
	JSL::gnuplot gp;
	double maxObs = 0;
	PlotModel(gp,Histogram,model,maxObs,true,"Baseline");
	
	Optimise(model,100);
	PlotModel(gp,Histogram,model,maxObs,false,"Fine-Tuned");
	gp.Show();
	return model;
}