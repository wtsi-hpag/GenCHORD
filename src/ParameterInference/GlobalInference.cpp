#include "GlobalInference.h"

int L_glob;
struct Helper
{
	int N;
	int kdim;
	std::vector<double> ws;

	std::vector<JSL::Vector> hs; //weight 
	JSL::Vector hSums;

	JSL::Vector es; //weight 
	Helper(const std::vector<int> & Nks, int Qmax)
	{
		kdim = Nks.size();

		hSums = JSL::Vector(kdim);
		hs = std::vector<JSL::Vector>(Qmax+1,JSL::Vector(kdim));
		ws.resize(Qmax+1);
		es = JSL::Vector(kdim);
		N = 0;
		for (int k = 0; k < kdim; ++k)
		{
			N+= Nks[k];
		}
	}

	void SetWeights()
	{
		double s= 0;
		for (int i = 0; i < ws.size(); ++i)
		{
			int ddip = abs(i-2);
			ws[i] = pow(10,-ddip);
			s+=ws[i];
		}
		for (int i = 0; i < ws.size(); ++i)
		{
			ws[i] /= s;
		}
	}
	void NudgeWeights(double nu)
	{
		double s= 0;
		for (int i = 0; i < ws.size(); ++i)
		{
			int ddip = abs(i-2);
			double t1= pow(10,-ddip);
			double mem = 0.5;
			ws[i] = mem * ws[i] + (1.0 - mem) * t1;

			if (i * nu > kdim*1.2) //kill off nodes which are outside the given space
			{
				ws[i] /= 1e6;
			}

			s+=ws[i];
		}
		for (int i = 0; i < ws.size(); ++i)
		{
			ws[i] /= s;
		}
	}
	void ComputeHs(ProbabilityModel & p)
	{
		hSums *=0;
		hSums -= 999999;
		for (int q = 0; q < ws.size(); ++q)
		{
			double logWq = log(ws[q]);
			double norm = -999999;
			
			for (int k = 0; k < kdim; ++k)
			{
				hs[q][k] = p.logP(k,q) + logWq;
				hSums[k] = ale(hSums[k],hs[q][k]);
			}
		}
		for (int q = 0; q < ws.size(); ++q)
		{
			
			hs[q] =hs[q] - hSums;
		}
	}

	void ComputeEtas(ProbabilityModel & p)
	{
		double sigSum = -99999;
		double noiseSum = -99999;
		double noiseWeight = log(p.NoiseWeight);
		double signalWeight = log(1.0 - p.NoiseWeight);
		
		double runSum = 0;
		for (int k =0; k < kdim; ++k)
		{
			double e = noiseWeight + p.getLogNoise(k);
			double s = e;
			for (int q = 0; q < hs.size(); ++q)
			{
				s= ale(s,signalWeight + log(ws[q]) + p.getLogSignal(k,q)) ;
			}
			es[k] = e - s;
		}
	}
};

void ParameterRelaxation(ProbabilityModel & p, const std::vector<int> & Nks, int kMax, int Qmax,int res, Helper & assist,bool nudge)
{
	// p.NoiseMean = 15;//p.SignalMean;
	// p.NoiseSigma =30;//2*p.NoiseMean;

	if (nudge)
	{
		double w1 = 0.2;
		double mod1 = w1 * p.NoiseMean + (1.0 - w1) * 10;
		double mod2 = w1 * p.NoiseSigma + (1.0 - w1) * 60;
		double mod3 = w1 * p.NoiseWeight + (1.0 - w1) * exp(-2);
		p.SetNoiseParameters(mod1,mod2,mod3);

		assist.NudgeWeights(p.SignalMean);
	}
	for (int l = 0; l < res; ++l)
	{

		p.SetGrids();
		// //determine the updated weights
		
		assist.ComputeHs(p);
		assist.ComputeEtas(p);
		double diff = 0;

		for (int q =0; q <=Qmax; ++q)
		{
			double newW = 0;
			for (int i = 0; i < assist.kdim; ++i)
			{
				newW += exp(assist.hs[q][i])/assist.N * Nks[i];
			}
			diff += pow((newW - assist.ws[q])/(assist.ws[q]),2);
			assist.ws[q] = newW;		
		}

		//determine updated noise weight
		
		double newGamma=0;
		double newNoiseMean = 0;
		double newNoiseSigma = 0;
		for (int k = 0; k < assist.kdim; ++k)
		{
			newGamma += Nks[k] * exp(assist.es[k])/assist.N;
		}
		for (int k = 0; k < assist.kdim; ++k)
		{
			double p = exp(assist.es[k]) * Nks[k]/(assist.N*newGamma);
			newNoiseMean +=  p*k;
			newNoiseSigma += p*k*k;
		}
		newNoiseSigma = sqrt(newNoiseSigma - newNoiseMean*newNoiseMean);

		double oldMu = p.NoiseMean;
		double oldSigma = p.NoiseSigma;
		double oldGamma =p.NoiseWeight;
		p.SetNoiseParametersFromObserved(newNoiseMean,newNoiseSigma,newGamma);
		diff += pow((oldMu-p.NoiseMean)/oldMu,2) + pow((oldSigma -p.NoiseSigma)/oldSigma,2) + pow((oldGamma - p.NoiseWeight)/oldGamma,2);


		diff = sqrt(diff/(3 + assist.ws.size()));

		L_glob = l;
		if (l > res/3 && diff < 0.02)
		{
			break;
		}	
	}
}



std::vector<double> probContainer;
double ComputeScore(ProbabilityModel & p, const std::vector<int> & Nks, std::vector<double> & ws, int kMax)
{
	p.SetGrids();
	if (probContainer.size() != Nks.size())
	{
		probContainer.resize(kMax+1);
	}
	p.GlobalLogPrediction(probContainer,ws);
	double score = 0;
	int nSum = 0;
	for (int k = 0; k <Nks.size(); ++k)
	{
		score += Nks[k] * probContainer[k];
		nSum += Nks[k];
	}
	double prior =0;

	double maxw = 0;
	int maxq = 0;
	for (int q = 0; q < ws.size(); ++q)
	{
		if (ws[q] > maxw)
		{
			maxq = q;
			maxw = ws[q];
		}
	}
	if (maxq != 2)
	{
		double nNonDiploid =  (ws[maxq] - ws[2]) * nSum;
		double nonDiploidPenalty = -0.1*nNonDiploid;
		prior += nonDiploidPenalty;
	}
	if (p.NoiseWeight > 0.25)
	{
		double nNoiseExcess = (p.NoiseWeight - 0.25) * nSum;
		double noisePenalty = -0.1*nNoiseExcess;
		prior += noisePenalty;
	}
	return score + prior;
}


void NormaliseModel(ProbabilityModel & p, const Data & data, const Settings & settings, int kMax,int Qmax)
{
	
	p.SetDimensionality(kMax,Qmax);
	p.ComputeNoiseConversionChart();
	std::vector<int> Nks(kMax+1,0);
	for (int c = 0; c < data.Chromosomes.size(); ++c)
	{
		for (int i = 0; i < data.Chromosomes[c].Counts.size(); ++i)
		{
			int k = data.Chromosomes[c].Counts[i];
			if (k<= kMax)
			{
				Nks[k] += 1;
			}
		}
	}
	Helper assist(Nks,Qmax);
	
	assist.SetWeights();

	p.NoiseMean = 9;
	p.NoiseSigma = 1;
	p.NoiseWeight = exp(-2);
	// p.SignalMean = 25;
	// p.SignalSigma = 2;
	

	// std::vector<double> mus = JSL::Vector::linspace(15,25,151);
	std::vector<double> mus = JSL::Vector::linspace(data.Mean/3,data.Mean/1.5,60);
	std::vector<double> sigmas = JSL::Vector::linspace(1,25,40);
	
	int N = JSL::Vector(Nks).Sum();
	double perfectScore = -N*log(N);
	for (int k = 0; k <= kMax; ++k)
	{
		if (Nks[k] > 0)
		{
			perfectScore += Nks[k] * log(Nks[k]);
		}
	}

	bool found = false;
	double bestScore = -9e10;
	double bestMu;
	double bestSig;
	double bestNoiseMean;
	double bestNoiseSig;
	double bestNoiseWeight;
	std::vector<double> bestWeights;
	double diff = 0;

	double shiftMem = 1.0;
	double mu_e = data.Mean;
	double sig_e = mu_e;
	double lMem = 0.99;
	double L = 0;

	std::vector<std::vector<double>> scores(sigmas.size(),std::vector<double>(mus.size(),0.0));
	JSL::ProgressBar<2> PB(mus.size(),sigmas.size());
	int n = 0;
	Log("\tDetermining global distribution parameters\n")
	PB.SetName(0,"\t");
	PB.SetName(1,"\t");
	for (int i = 0; i < mus.size(); ++i)
	{
		// p.NoiseMean = 10;// mus[i];
		// p.NoiseSigma = 1;
		for (int j = 0; j < sigmas.size(); ++j)
		{
			
			p.SignalMean = mus[i];
			p.SignalSigma = sigmas[j];
			ParameterRelaxation(p,Nks,kMax,Qmax,100,assist,true);
			double score = ComputeScore(p,Nks,assist.ws,kMax) - perfectScore;
			scores[j][i] = abs(score);
			if (!found || score > bestScore)
			{
				found = true;
				bestScore = score;
				bestMu = p.SignalMean;
				bestSig = p.SignalSigma;
				bestNoiseMean = p.NoiseMean;
				bestNoiseSig = p.NoiseSigma;
				bestNoiseWeight = p.NoiseWeight;
				bestWeights = assist.ws;
			}
			PB.Update(i,j);
			PB.SetName(1,"\t"+std::to_string(int(bestScore)));
		}
	}

	
	p.SetSignalParameters(bestMu,bestSig);
	p.SetNoiseParameters(bestNoiseMean,bestNoiseSig,bestNoiseWeight);
	assist.ws = bestWeights;
	ParameterRelaxation(p,Nks,kMax,Qmax,250,assist,false);
	p.SetGrids();
	double trueScore = ComputeScore(p,Nks,assist.ws,kMax) - perfectScore;
	Log("\tOptimal signal model has nu=" << bestMu << ", sigma=" << bestSig << " with a fraction " << assist.ws[2] *100 << "% of the data assigned as diploid\n");
	Log("\tOptimal noise model has mu=" << p.NoiseMean << ", sigma=" << p.NoiseSigma << " and has weight " << p.NoiseWeight << "\n");
	Log("\tOptimal combined model as a deviation score of " << trueScore << "\n");

	
	JSL::gnuplot gp;
	std::vector<int> ks = JSL::Vector::intspace(0,kMax,1);

	gp.SetMultiplot(2,1);
	namespace lp = JSL::LineProperties;
	std::vector<double> probs(ks.size());
	p.GlobalLogPrediction(probs,assist.ws);
	
	double logN = log(assist.N);
	for (int k = 0; k < probs.size(); ++k)
	{
		probs[k] = exp(logN + probs[k]);
		if (probs[k] < 0.1)
		{
			probs[k] = 0;
		}
	}
	gp.Plot(ks,Nks,lp::Legend("Data"));
	gp.Plot(ks,probs,lp::Legend("Probability Distribution"),lp::PenSize(5));
	gp.WindowSize(2000,1400);
	gp.SetTerminal("pngcairo");
	gp.SetOutput(settings.OutputDirectory + "/probability_model.png");
	gp.SetXLabel("k");
	gp.SetYLabel("p_k");
	gp.SetTitle("Inferred Probability Model");
	gp.SetAxis(1);
	gp.WindowSize(2000,1400);
	gp.SetFontSize(JSL::Fonts::Global,25);
	gp.SetXLabel("k");
	gp.SetYLabel("p_k (log scale)");
	
	gp.Plot(ks,Nks,lp::Legend("Data"));
	gp.Plot(ks,probs,lp::Legend("Probability Distribution"),lp::PenSize(5));
	gp.SetYLog(true);

	std::vector<double> ps(ks.size());
	for (int k  =0; k < ks.size(); ++k)
	{
		ps[k] = exp(p.logNoiseWeight + p.getLogNoise(k) + logN);
		if (ps[k] < 0.1)
		{
			ps[k] = 0;
		}
	}
	gp.SetAxis(0);
	gp.SetLegend(true);
	gp.Plot(ks,ps,JSL::LineProperties::PenType(JSL::Dash),JSL::LineProperties::PenSize(3),JSL::LineProperties::Legend("Error Component"));
	gp.SetAxis(1);
	gp.Plot(ks,ps,JSL::LineProperties::PenType(JSL::Dash),JSL::LineProperties::PenSize(2));
	for (int q = 0; q < assist.ws.size(); ++q)
	{
		for (int k  =0; k < ks.size(); ++k)
		{
			ps[k] = exp(p.logSignalWeight + p.getLogSignal(k,q) + logN + log(assist.ws[q]));
			if (ps[k] < 0.1)
			{
				ps[k] = 0;
			}
		}
		gp.SetAxis(0);
		gp.Plot(ks,ps,JSL::LineProperties::PenType(JSL::Dotted),JSL::LineProperties::PenSize(2));
		gp.SetAxis(1);
		gp.Plot(ks,ps,JSL::LineProperties::PenType(JSL::Dotted),JSL::LineProperties::PenSize(2));
	}

	gp.Show();

	JSL::gnuplot gp2;
	gp2.WindowSize(2000,1400);
	gp2.SetFontSize(JSL::Fonts::Global,25);
	gp2.Map(mus, sigmas,scores);
	gp2.SetCBLog(true);
	gp2.SetTerminal("pngcairo");
	gp2.SetXLabel("mu");
	gp2.SetYLabel("sigma");
	gp2.SetOutput(settings.OutputDirectory + "/model_parameters.png");
	gp2.Show();
}
