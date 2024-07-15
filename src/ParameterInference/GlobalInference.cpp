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
	std::vector<std::vector<double>> mean_posterior;
	std::vector<double> noise_zs;
	std::vector<double> noise_ms;
	std::vector<double> noise_vs;
	Helper(const std::vector<int> & Nks, int Qmax,int noiseRes)
	{
		kdim = Nks.size();

		signal_posterior = std::vector<std::vector<double>>(Qmax+1,JSL::Vector(kdim));
		mean_posterior = signal_posterior;
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
			double ddip = sqrt(abs(i-2));
			ws[i] = pow(2,-ddip);
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
				signal_posterior[q][k] = p.logP(k,q) + logWq - log_pk;
				mean_posterior[q][k] = p.getLogSignal(k,q) + logWq - log_sk;
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
			double prior_k = k;
			double crit = 2*p.SignalMean;
			if (k > crit)
			{
				prior_k *= k * 10.0/(3*crit);
			}
			for (int q= 0; q < ws.size(); ++q)
			{
				double d = (k - q * p.SignalMean)/(p.SignalSigma*0.3);
				prior_k += 0.00001*exp(-0.5*d*d);
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
			NoiseGradient[r] = exp(temp[r]) - exp(p.NoiseComponents[r] + globalSum) - 1e-4*N*priorTerm;
		}

	}

	void UpdateNoiseComponents(ProbabilityModel & p,int l)
	{
		double alpha = 0.1 * pow(0.98,l);

		double b1 = 0.7;
		double b2 = 0.9;
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

	double gammaCap = 1e-4;
	double gammaGain = 3;
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
			double meanW = 0;
			for (int i = 0; i < assist.kdim; ++i)
			{
				double p = exp(assist.signal_posterior[q][i])/assist.N * Nks[i];
				newW += p;
				double pM = exp(assist.mean_posterior[q][i])/assist.N * Nks[i];
				meanW += p;
				meanMod += pM*i; 

			}
			wSum += newW;
			meanMod /= meanW;
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
			double gammaMem = 0.;
			newGamma = gammaMem * exp(p.logNoiseWeight) + (1.0 - gammaMem)*newGamma;
			newGamma = std::min(gammaCap,newGamma);
			if (gammaCap < 0.07)
			{
				gammaCap *= gammaGain;
			}
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


std::vector<double> NormaliseSet(ProbabilityModel & p, std::vector<int> Nks,int kMax, int Qmax,std::vector<double> & mus, std::vector<double> & sigmas)
{
	Helper assist(Nks,Qmax,p.NoiseResolution);
	assist.SetWeights(p);

	double perfectScore = -assist.N*log(assist.N);
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
	JSL::ProgressBar<2> PB(mus.size(),sigmas.size());
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
			ParameterRelaxation(p,Nks,kMax,Qmax,75,assist,true);
			double score = ComputeScore(p,Nks,assist.ws,kMax) - perfectScore;
			if (!found || score > bestScore)
			{
				found = true;
				bestScore = score;
				bestMu = p.SignalMean;
				bestSig = p.SignalSigma;
			}
			PB.Update(i,j);
			PB.SetName(1,"\t"+std::to_string(int(bestScore)));
		}
	}

	p.SetSignalParameters(bestMu,bestSig);
	ParameterRelaxation(p,Nks,kMax,Qmax,250,assist,false);
	return assist.ws;
}

void PlotDistribution(ProbabilityModel & p, std::vector<double> ws,const std::vector<int> & Nks, const Settings & settings,std::string subname)
{
	JSL::gnuplot gp;
	std::vector<int> ks = JSL::Vector::intspace(0,p.MaxK,1);

	gp.SetMultiplot(2,1);
	namespace lp = JSL::LineProperties;
	std::vector<double> probs(ks.size());
	p.GlobalLogPrediction(probs,ws);
	
	double logN = -99999999;
	for (int k =0; k < Nks.size(); ++k)
	{
		logN = ale(logN,log(Nks[k]));
	}
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
	gp.WindowSize(1000,700);
	gp.SetTerminal("pngcairo");
	gp.SetOutput(settings.OutputDirectory + "/Distributions/" + subname + "_probability_model.png");
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
	for (int q = 0; q < ws.size(); ++q)
	{
		for (int k  =0; k < ks.size(); ++k)
		{
			ps[k] = exp(p.logSignalWeight + p.getLogSignal(k,q) + logN + log(ws[q]));
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
}

std::vector<ProbabilityModel> NormaliseModel(ProbabilityModel & p, const Data & data, const Settings & settings, int kMax,int Qmax)
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
	JSL::mkdir(settings.OutputDirectory + "/Distributions/");
	Log("\tDetermining global distribution parameters\n")
	std::vector<double> mus_fullScan = JSL::Vector::linspace(0.4,1,101);
	Log("\tScanning between " << mus_fullScan[0] << "  " << mus_fullScan[mus_fullScan.size()-1] << std::endl;);
	std::vector<double> sigmas_fullScan = JSL::Vector::linspace(0.001,1,8);
	
	auto ws = NormaliseSet(p,Nks,kMax,Qmax,mus_fullScan,sigmas_fullScan);
	PlotDistribution(p,ws,Nks,settings,"global");
	auto pGlobal = p;
	std::vector<double> mus_SubScan = JSL::Vector::linspace(p.SignalMean*0.5,p.SignalMean*1.5,21);
	std::vector<double> sigmas_SubScan = JSL::Vector::linspace(p.SignalSigma*0.5,p.SignalSigma*1.5,21);
	std::vector<ProbabilityModel> Distributions;
	for (int c =0; c < data.Chromosomes.size(); ++c)
	{
		Log("\tDetermining parameters for chromosome " << data.Chromosomes[c].Name << std::endl;)

		std::fill(Nks.begin(),Nks.end(),0);
		for (int i = 0; i < data.Chromosomes[c].Counts.size(); ++i)
		{
			int k = data.Chromosomes[c].Counts[i];
			if (k<= kMax)
			{
				Nks[k] += 1;
			}
		}

		ws = NormaliseSet(p,Nks,kMax,Qmax,mus_SubScan,sigmas_SubScan);
		PlotDistribution(p,ws,Nks,settings,"sub"+data.Chromosomes[c].Name);
		p.SetGrids();
	
		Distributions.push_back(p);
	}

	
	

	return Distributions;
}
