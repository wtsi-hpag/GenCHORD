#include "GlobalInference.h"

struct Helper
{
	int N;
	int kdim;
	std::vector<double> ws;
	int nr;
	std::vector<std::vector<double>> signal_posterior; //weight 

	std::vector<double> NoiseGradient;
	std::vector<double> noise_posterior; //noise weight helper

	std::vector<double> noise_zs;
	std::vector<double> noise_ms;
	std::vector<double> noise_vs;
	Helper(const std::vector<int> & Nks, int Qmax,int noiseRes)
	{
		kdim = Nks.size();

		signal_posterior = std::vector<std::vector<double>>(Qmax+1,JSL::Vector(kdim));
		noise_posterior = std::vector<double>(kdim);
		ws.resize(Qmax+1);

		N = 0;
		nr = noiseRes;
		NoiseGradient.resize(nr);
		noise_zs.resize(nr);
		noise_ms.resize(nr);
		noise_vs.resize(nr);

		for (int k = 0; k < kdim; ++k)
		{
			N+= Nks[k];
		}
	}

	void SetWeights(ProbabilityModel & p)
	{
		std::fill(noise_ms.begin(),noise_ms.end(),0.0);
		std::fill(noise_vs.begin(),noise_vs.end(),0.0);
		double s= 0;
		// for (int r = 0; r < p.NoiseComponents.size(); ++r)
		// {
		// 	s+= exp(p.NoiseComponents[r]);
		// }
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

	void ComputePosteriors(ProbabilityModel & p)
	{
		for (int k = 0; k < kdim; ++k)
		{
			double log_pk = -999999999;
			double log_sk = -9999999999;
			for (int q = 0; q < ws.size(); ++q)
			{
				double signal = log(ws[q]) + p.getLogSignal(k,q);
				double total = ale(p.logSignalWeight + signal,p.logNoiseWeight+p.getLogNoise(k));
				log_pk= ale(log_pk,log(ws[q])+p.logP(k,q));
				log_sk = ale(log_sk,signal);
				// std::cout << "\t" <<k << " " << q << " " << p.logP(k,q) << "  " << log_pk << std::endl;
			}

			if (std::isnan(log_pk) || std::isinf(log_pk))
			{
				std::cout << p.logSignalWeight << std::endl;
				std::cout << JSL::Vector(p.NoiseComponents) << std::endl;
				std::cout << "SHIT" << std::endl;
				exit(5);
			}

			for (int q = 0; q < ws.size(); ++q)
			{
				double logWq = log(ws[q]);
				signal_posterior[q][k] = p.getLogSignal(k,q) + logWq - log_sk;
			}
			// for (int r = 0; r < noise_posterior.size(); ++r)
			// {
			noise_posterior[k] = p.logNoiseWeight + p.getLogNoise(k) - log_pk;
			// }
		}
		// hSums *=0;
		// hSums -= 999999;
		// for (int q = 0; q < ws.size(); ++q)
		// {
		// 	double logWq = log(ws[q]);
		// 	double norm = -999999;
			
		// 	for (int k = 0; k < kdim; ++k)
		// 	{
		// 		hs[q][k] = p.logP(k,q) + logWq;
		// 		hSums[k] = ale(hSums[k],hs[q][k]);
		// 	}
		// }
		// for (int q = 0; q < ws.size(); ++q)
		// {
			
		// 	hs[q] =hs[q] - hSums;
		// }
	}

	void ComputeNoiseGradients(ProbabilityModel & p,const std::vector<int> & Nks)
	{
		std::vector<double> temp(nr,-99999999);
		std::vector<double> prior(nr,0);
		double globalSum = -9999999999;
		double priorSum = 0;
		for (int k =0; k < kdim; ++k)
		{
			double p_k =-99999999;
			double prior_k = 0;
			for (int q= 0; q < ws.size(); ++q)
			{
				double d = (k - q * p.SignalMean)/(p.SignalSigma*0.5);
				prior_k += exp(-0.5*d*d);
				p_k = ale(p_k,log(ws[q]) + p.logP(k,q));
			}
			int r = p.RFromK(k);
			double term = log(Nks[k]) - p_k - p.logNoiseWidth[r] + p.NoiseComponents[r];
			double priorTerm = prior_k * exp(p.NoiseComponents[r]);
			globalSum = ale(globalSum,term );
			temp[r] = ale(temp[r],term);
			prior[r] += priorTerm;
			priorSum += priorTerm;
		}

		for (int r =0 ; r < nr; ++r)
		{
			double priorTerm = prior[r] - exp(p.NoiseComponents[r]) * priorSum;
			NoiseGradient[r] = exp(temp[r]) - exp(p.NoiseComponents[r] + globalSum) - 4e-2*N*priorTerm;
		}

	}

	void UpdateNoiseComponents(ProbabilityModel & p,int l)
	{
		double alpha = 0.1;
		double b1 = 0.7;
		double b2 = 0.99;
		double norm = -999999999;
		double c1 = 1.0/(1.0 - pow(b1,l+1));
			double c2 = 1.0/(1.0 - pow(b2,l+1));
		for (int r = 0; r < nr; ++r)
		{	
			noise_zs[r] = p.NoiseComponents[r];

			noise_ms[r] = b1*noise_ms[r] + (1.0 - b1)*NoiseGradient[r];
			noise_vs[r] = b2*noise_vs[r] + (1.0 - b2) * NoiseGradient[r]*NoiseGradient[r];
			
			double step = alpha * noise_ms[r]*c1/sqrt(noise_vs[r]*c2 + 1e-20 );
			noise_zs[r] += step;
			norm = ale(norm,noise_zs[r]);
		}

		for (int r =0; r < nr; ++r)
		{
			p.NoiseComponents[r] = noise_zs[r] - norm;
		}

	}
};

void ParameterRelaxation(ProbabilityModel & p, const std::vector<int> & Nks, int kMax, int Qmax,int res, Helper & assist,bool nudge)
{
	p.ClearContamination();
	p.ResetNoise();
	assist.SetWeights(p);
	std::vector<double> newNoiseParams(p.NoiseResolution);
	double noiseCap = 1e-1;
	int regainSteps = (res/2 - res/5);
	double noiseGain = exp(-1.0/regainSteps * log(noiseCap));
	int timeSinceKicked =0;
	double alpha = 0.1;

	double gammaCap = 1e-2;
	double gammaGain = 1.5;
	for (int l = 0; l < res; ++l)
	{

		p.SetGrids();
		// //determine the updated weights
		// std::cout << "posterior" << std::endl;
		assist.ComputePosteriors(p);
		
		double diff = 0;
		double wSum = 0;
		for (int q =0; q <=Qmax; ++q)
		{
			double newW = 0;
			double meanMod = 0;
			double newGamma = 0;
			for (int i = 0; i < assist.kdim; ++i)
			{
				double p = exp(assist.signal_posterior[q][i])/assist.N * Nks[i];
				meanMod += p*i; 
				newW += p;

			}
			wSum += newW;
			meanMod /= newW;
			if (q!=2)
			{
				p.SetContamination(q,meanMod);
			}
			diff += pow((newW - assist.ws[q])/(assist.ws[q]),2);
			assist.ws[q] = newW;
		}
		
		// std::cout << newGamma << std::endl;

		int lWait = 0;
		if (l >= lWait)
		{
			double newGamma = 0;
			double posSum = 0;
			for (int i = 0; i < assist.kdim; ++i)
			{
				newGamma += exp(assist.noise_posterior[i])/assist.N * Nks[i];
			}
			double gammaMem = 0.7;
			newGamma = gammaMem * exp(p.logNoiseWeight) + (1.0 - gammaMem)*newGamma;
			newGamma = std::min(gammaCap,newGamma);
			gammaCap *= gammaGain;
			p.logNoiseWeight = log(newGamma);
			p.logSignalWeight = log(1.0 - newGamma);

			assist.ComputeNoiseGradients(p,Nks);
			assist.UpdateNoiseComponents(p,l-lWait);
		}
		// if (l > 0)
		// {
		// 	double gamma = 0;
		// 	for (int r =0; r < p.NoiseResolution; ++r)
		// 	{
		// 		double newW = 0;
		// 		for (int i = 0; i < assist.kdim; ++i)
		// 		{
		// 			newW +=exp(assist.noise_posterior[r][i])/assist.N * Nks[i];
		// 		}
		// 		newNoiseParams[r] = newW;
		// 		gamma += newW;
				
		// 	}
			
		// 	double norm = gamma + wSum;
		// 	//enforce normalisation
		// 	for (int q = 0; q < assist.ws[q]; ++q)
		// 	{
		// 		assist.ws[q]/=norm;
		// 	}


		// 	p.SetNoiseComponents(newNoiseParams);
		// }

		// double wSum = JSL::Vector(assist.ws).Sum();
		// for (int i = 0; i < p.NoiseResolution; ++i)
		// {
		// 	wSum += 
		// }

		// exit(5);
		// 	// if (gamma > 0.1)
		// 	// {
		// 	// 	double corrector = log(0.1) - log(gamma);
		// 	// 	for (int )
		// 	// }
		// }

		// diff = sqrt(diff/(3 + assist.ws.size()));

		// if (l > res/3 && diff < 0.02)
		// {
		// 	break;
		// }	
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
	// double noiseCrit = 0.1;
	// if (p.NoiseWeight > noiseCrit)
	// {
	// 	double nNoiseExcess = (p.NoiseWeight - noiseCrit) * nSum;
	// 	double noisePenalty = -0.1*nNoiseExcess;
	// 	prior += noisePenalty;
	// }
	return score + prior;
}


void NormaliseModel(ProbabilityModel & p, const Data & data, const Settings & settings, int kMax,int Qmax)
{
	p.SetDimensionality(kMax,Qmax);
	
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
	Helper assist(Nks,Qmax,p.NoiseResolution);
	
	assist.SetWeights(p);

	
	// std::vector<double> mus = JSL::Vector::linspace(45,65,151);
	std::vector<double> mus = JSL::Vector::linspace(data.Mean/2.5,data.Mean/1.5,50);
	std::vector<double> sigmas = JSL::Vector::linspace(1,20,20);
	
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
			ParameterRelaxation(p,Nks,kMax,Qmax,50,assist,true);
			double score = ComputeScore(p,Nks,assist.ws,kMax) - perfectScore;
			scores[j][i] = abs(score);
			if (!found || score > bestScore)
			{
				found = true;
				bestScore = score;
				bestMu = p.SignalMean;
				bestSig = p.SignalSigma;
				bestWeights = assist.ws;
			}
			PB.Update(i,j);
			PB.SetName(1,"\t"+std::to_string(int(bestScore)));
		}
	}
	p.ClearContamination();
	p.SetSignalParameters(bestMu,bestSig);
	assist.ws = bestWeights;
	ParameterRelaxation(p,Nks,kMax,Qmax,250,assist,false);
	p.SetGrids();
	double trueScore = ComputeScore(p,Nks,assist.ws,kMax) - perfectScore;
	Log("\tOptimal signal model has nu=" << bestMu << ", sigma=" << bestSig << " with a fraction " << assist.ws[2] *100 << "% of the data assigned as diploid\n");
	Log("\tOptimal noise model has a weight " << exp(p.logNoiseWeight)<< "\n");
	Log("\tOptimal combined model as a deviation score of " << trueScore << "\n");

	
	JSL::gnuplot gp;
	std::vector<int> ks = JSL::Vector::intspace(0,kMax,1);

	gp.SetMultiplot(2,1);
	namespace lp = JSL::LineProperties;
	std::vector<double> probs(ks.size());
	p.GlobalLogPrediction(probs,assist.ws);
	
	double logN = log(assist.N);
	double s =-999999;
	for (int k = 0; k < probs.size(); ++k)
	{
		s= ale(probs[k],s);
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
		ps[k] = exp(p.getLogNoise(k) + logN + p.logNoiseWeight);
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
