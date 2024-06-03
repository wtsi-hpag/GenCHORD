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
			ws[i] = pow(3,-ddip);
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
				hs[q][k] = p.logP(k,q,true) + logWq;
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
			double e = noiseWeight + p.logNoise(k) - p.NoiseNormalisation;
			double s = -999999;
			for (int q = 0; q < hs.size(); ++q)
			{
				s= ale(s,signalWeight + log(ws[q]) + p.logSignal(k,q) - p.Normalisation[q]) ;
			}
			es[k] = e - ale(e,s);
		}
	}
};

std::vector<double> ParameterRelaxation(ProbabilityModel & p, const std::vector<int> & Nks, int kMax, int Qmax,int res)
{
	// p.NoiseMean = 15;//p.SignalMean;
	// p.NoiseSigma =30;//2*p.NoiseMean;

	// std::vector<double> weights(Qmax+1,1);
	// for (int i = 0; i <= Qmax; ++i)
	// {
	// 	double ddip = abs(i-2);
	// 	weights[i] = weights[2]/pow(5,ddip);
	// }
	p.NoiseWeight = exp(-2);
	Helper assist(Nks,Qmax);
	assist.SetWeights();
	for (int l = 0; l < res; ++l)
	{

		p.Normalise(kMax,assist.ws);
		//determine the updated weights
		
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
			double p = exp(assist.es[k]) * Nks[k]/(assist.N * newGamma);
			newNoiseMean +=  p*k;
			newNoiseSigma += p*k*k;
		}
		newNoiseSigma = sqrt(newNoiseSigma - newNoiseMean*newNoiseMean);

		double oldMu = p.NoiseMean;
		double oldSigma = p.NoiseSigma;
		double oldGamma =p.NoiseWeight;
		p.SetNoiseParameters_Observed(newGamma,newNoiseMean,newNoiseSigma);

		diff += pow((oldMu-p.NoiseMean)/oldMu,2) + pow((oldSigma -p.NoiseSigma)/oldSigma,2) + pow((oldGamma - p.NoiseWeight)/oldGamma,2);


		diff = sqrt(diff);

		L_glob = l;
		if (l > res/5 && diff < 0.03)
		{
			break;
		}
		
		// std::cout << "\t" << l << " " << newGamma << "  " << p.NoiseMean << "  " << p.NoiseSigma << std::endl;
	}
	return assist.ws;
}



double ComputeScore(ProbabilityModel & p, const std::vector<int> & Nks, std::vector<double> & ws, int kMax,bool usePriors)
{
	p.Normalise(kMax,ws);
	auto probs = p.GlobalPrediction(kMax,ws);
	double score = 0;
	int nSum = 0;
	for (int k = 0; k <Nks.size(); ++k)
	{
		score += Nks[k] * log(std::max(1e-100,probs[k]));
		nSum += Nks[k];
	}
	double prior =0;
	double nNonDiploid =  (1.0 - ws[2]) * nSum;
	double nonDiploidPenalty = -0.00*nNonDiploid;

	prior += nonDiploidPenalty;
	// if (p.NoiseWeight > 0.25)
	// {
	// 	double nNoise = (p.NoiseWeight-0.25) * nSum;
	// 	double noisePenalty = -0.1 * nNoise;
	// 	prior += noisePenalty;
	// }

	// std::cout << p.NoiseMean << "  " << p.NoiseSigma << "  " << score << std::endl;
	return score + prior;
}


void GlobalInference(ProbabilityModel & p, const Data & data, const Settings & settings, int kMax,int Qmax)
{
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
	p.ComputeNoiseConversionChart(kMax);
	int ModelDimension = p.Dimensionality;
	int TotalDimensionality = ModelDimension + Qmax+1;
	
	JSL::Vector Position(TotalDimensionality);
	
	p.SignalMean = 25;
	p.SignalSigma = 2;
	

	// std::vector<double> mus = JSL::Vector::linspace(data.Mean/3,data.Mean,121);
	std::vector<double> mus = JSL::Vector::linspace(12,28,35);
	std::vector<double> sigmas = JSL::Vector::linspace(20,43,48);
	std::cout << "Testing " << data.Mean/3 << "  < mu < " << data.Mean * 1.1 << std::endl;
	
	int N = JSL::Vector(Nks).Sum();
	double perfectScore = -N*log(N);
	for (int k = 0; k <= kMax; ++k)
	{
		if (Nks[k] > 0)
		{
			perfectScore += Nks[k] * log(Nks[k]);
		}
	}

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
			p.NoiseMean = 5;
			p.NoiseSigma = 3;
			p.SignalMean = mus[i];
			p.SignalSigma = sigmas[j];
			auto ws = ParameterRelaxation(p,Nks,kMax,Qmax,150);
			L = lMem * L + (1.0 - lMem) * L_glob;
			double score = ComputeScore(p,Nks,ws,kMax,true) - perfectScore;
			scores[j][i] = abs(score);
			// std::cout << p.SignalMean << " " << p.SignalSigma << "  " << score << std::endl;
			if (score > bestScore)
			{
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
	p.NoiseMean = 5;
	p.NoiseSigma = 3;
	auto ws = ParameterRelaxation(p,Nks,kMax,Qmax,150);
	double trueScore = ComputeScore(p,Nks,ws,kMax,false);
	std::cout << "Inferred parameters " << std::endl;
	std::cout << "\tmu_e: \t" << p.NoiseMean << "\n\tsig_e:\t" << p.NoiseSigma << "\n\tgamma:\t" << p.NoiseWeight << std::endl;
	std::cout << "\tws: \t" << JSL::Vector(ws) << std::endl;
	std::cout << "\tInferred score " << ComputeScore(p,Nks,ws,kMax,true) - perfectScore<< std::endl;
	JSL::gnuplot gp;
	std::vector<int> ks;
	ks = JSL::Vector::intspace(0,kMax,1);

	gp.SetMultiplot(2,1);
	gp.Scatter(ks,Nks);
	auto probs = p.GlobalPrediction(kMax,ws);
	
	for (int k = 0; k < probs.size(); ++k)
	{
		probs[k] *= N;
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
	std::ostringstream os;
	// os << "Old-Deforest (Gaussian) -- Statistical Deviation Score " << abs(trueScore - perfectScore);
	os << "New-Deforest (Negative Binomial) -- Statistical Deviation Score " << abs(trueScore - perfectScore);
	std::string title = os.str();
	std::cout << "TIT" << title << std::endl;
	gp.SetSuperTitle(title);
	gp.SetPersistence(true);
	gp.Show();

	JSL::gnuplot gp2;
	gp2.Map(mus, sigmas,scores);
	gp2.SetCBLog(true);
	gp2.Show();


	
	// auto prediction = p.GlobalPrediction(ks,assist.ws);
	// prediction =  JSL::Vector(prediction)/JSL::Vector(prediction).Sum();
	// std::cout << JSL::Vector(prediction).Sum() << std::endl;


	
	
	// JSL::gnuplot gp;
	// gp.SetMultiplot(2,1);
	// gp.Scatter(ks,Nks);
	// gp.Plot(ks,prediction);
	// // gp.SetXRange(0,100);
	// gp.SetYLog(true);

	// gp.SetAxis(1);
	// gp.Plot(ks,Nks);
	// gp.Plot(ks,prediction);
	// // gp.SetXRange(0,100);
	// gp.Show();
	


}
