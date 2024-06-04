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
			ws[i] = pow(4,-ddip);
			s+=ws[i];
		}
		for (int i = 0; i < ws.size(); ++i)
		{
			ws[i] /= s;
		}
	}
	void NudgeWeights()
	{
		double s= 0;
		for (int i = 0; i < ws.size(); ++i)
		{
			int ddip = abs(i-2);
			double t1= pow(4,-ddip);
			double mem = 0.7;
			ws[i] = mem * ws[i] + (1.0 - mem) * t1;
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

std::vector<double> ParameterRelaxation(ProbabilityModel & p, const std::vector<int> & Nks, int kMax, int Qmax,int res, Helper & assist)
{
	// p.NoiseMean = 15;//p.SignalMean;
	// p.NoiseSigma =30;//2*p.NoiseMean;

	double w1 = 0.6;
	double mod1 = w1 * p.NoiseMean + (1.0 - w1) * 9;
	double mod2 = w1 * p.NoiseSigma + (1.0 - w1) * 1;
	double mod3 = w1 * p.NoiseWeight + (1.0 - w1) * exp(-1);
	p.SetNoiseParameters(mod1,mod2,mod3);

	assist.NudgeWeights();
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
		if (l > res/5 && diff < 0.03)
		{
			break;
		}	
	}
	return assist.ws;
}



std::vector<double> probContainer;
double ComputeScore(ProbabilityModel & p, const std::vector<int> & Nks, std::vector<double> & ws, int kMax)
{
	p.SetGrids();
	probContainer.resize(kMax+1);
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
		double nonDiploidPenalty = -0.01*nNonDiploid;
		prior += nonDiploidPenalty;
	}
	return score + prior;
}


void GlobalInference(ProbabilityModel & p, const Data & data, const Settings & settings, int kMax,int Qmax)
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

	
	// p.SignalMean = 25;
	// p.SignalSigma = 2;
	

	// std::vector<double> mus = JSL::Vector::linspace(15,25,151);
	std::vector<double> mus = JSL::Vector::linspace(data.Mean/3,data.Mean,155);
	std::vector<double> sigmas = JSL::Vector::linspace(0.2,40,80);
	std::cout << "Testing " << mus[0] << "  < mu < " << mus[mus.size()-1] << std::endl;
	
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
	double diff = 0;

	double shiftMem = 1.0;
	double mu_e = data.Mean;
	double sig_e = mu_e;
	double lMem = 0.99;
	double L = 0;
	std::vector<std::vector<double>> scores(sigmas.size(),std::vector<double>(mus.size(),0.0));
	JSL::ProgressBar<2> PB(mus.size(),sigmas.size());
	int n = 0;
	for (int i = 0; i < mus.size(); ++i)
	{
		
		for (int j = 0; j < sigmas.size(); ++j)
		{
			p.NoiseMean = 1;// mus[i];
			p.NoiseSigma = 3;
			p.SignalMean = mus[i];
			p.SignalSigma = sigmas[j];
			auto ws = ParameterRelaxation(p,Nks,kMax,Qmax,18,assist);
			L = lMem * L + (1.0 - lMem) * L_glob;
			double score = ComputeScore(p,Nks,ws,kMax) - perfectScore;
			scores[j][i] = abs(score);
			// std::cout << p.SignalMean << " " << p.SignalSigma << "  " << score << std::endl;
			if (!found || score > bestScore)
			{
				found = true;
				// std::cout << "new best score " << score << std::endl;
				bestScore = score;
				bestMu = p.SignalMean;
				bestSig = p.SignalSigma;
			}
			PB.Update(i,j);
			++n;
			double l_Disp = L/(1.0 - pow(lMem,n));
			PB.SetName({std::to_string((int)bestScore),std::to_string(l_Disp)});
		}
	}

	std::cout << "Completed testing" << std::endl;
	std::cout << "\tmu = " << bestMu << "\n\tsigma=" << bestSig << std::endl;
	std::cout << "\tScore = " << bestScore << std::endl;
	
	p.SignalMean = bestMu;
	p.SignalSigma = bestSig;
	auto ws = ParameterRelaxation(p,Nks,kMax,Qmax,150,assist);
	p.SetGrids();
	double trueScore = ComputeScore(p,Nks,ws,kMax);
	std::cout << "Inferred parameters " << std::endl;
	std::cout << "\tmu_e: \t" << p.NoiseMean << "\n\tsig_e:\t" << p.NoiseSigma << "\n\tgamma:\t" << p.NoiseWeight << std::endl;
	std::cout << "\tws: \t" << JSL::Vector(ws) << std::endl;
	std::cout << "\tInferred score " << ComputeScore(p,Nks,ws,kMax) - perfectScore<< std::endl;
	
	
	JSL::gnuplot gp;
	std::vector<int> ks = JSL::Vector::intspace(0,kMax,1);

	gp.SetMultiplot(2,1);
	gp.Scatter(ks,Nks);
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
	gp.Plot(ks,probs);

	gp.SetAxis(1);
	gp.Plot(ks,Nks);
	gp.Plot(ks,probs);
	gp.SetYLog(true);
	gp.Show();

	JSL::gnuplot gp2;
	gp2.Map(mus, sigmas,scores);
	gp2.SetCBLog(true);
	gp2.Show();
}
