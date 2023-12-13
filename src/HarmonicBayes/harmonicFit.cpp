#include "harmonicFit.h"



Transitions AssignQs(double nu, Data & d, ErroredBinomial & EB,double alpha, int L, int accelerator,Distributor & dist)
{
	double score = 0;
	EB.SetBracket(nu);
	Transitions out(d.Chromosomes.size());
	out.Score = 0;
	dist.UpdateAssigner(&out,alpha,L,accelerator,nu);
	
	dist.Signal(2);
	for (int i = 0; i < dist.MainChromAssigment.size(); ++i)
	{
		int chrom = dist.MainChromAssigment[i];
		ChromosomeAssign(&out,chrom,d,alpha,&EB,L,accelerator,nu);
	}
	dist.Gather();
	// out.Score = score;
	// std::cout <<"Assigned" << std::endl;
	return out;
}

void HarmonicFit(Data & d, const Settings & settings)
{
	Log("Beginning search for harmonic parameters\n")
	double sigmaMin = settings.sigmaMin;
	double sigmaMax = settings.sigmaMax;
	int sigmaResolution = settings.sigmaResolution;
	double nuMin = floor(d.Mean/4);
	double nuMax = ceil(d.Mean);
	int nuResolution = settings.nuResolution;

	double deltaSig = (sigmaMax - sigmaMin)/(sigmaResolution -1);
	double deltaNu = (nuMax - nuMin)/(nuResolution - 1);

	double gammaMin = settings.gammaMin;
	double gammaMax = settings.gammaMax;
	int gammaResolution = settings.gammaResolution;
	double deltaGamma = (gammaMax - gammaMin)/(gammaResolution - 1);
	if (settings.Gamma > -1)
	{
		gammaMin = settings.Gamma;
		gammaMax = settings.Gamma;
		gammaResolution = 1;
		deltaGamma = 0;
	}
	
	logPloidyPrior = log(settings.PloidyPrior);
	Ploidy = settings.Ploidy;
	int accelerator = settings.Accelerator;
	int bounder = std::max(100.0,std::ceil(5*gammaMax));
	int qMax = settings.Qmax;
	int L = settings.L;
	int res = 0.25*qMax * nuMax;
	
	double alpha = settings.ContinuityPrior;

	Transitions best;
	
	best.Score = -999999999999999999;
	double worstScore = 0;

	Distributor dist(settings.ParallelWorkers,d);

	

	JSL::ProgressBar<3> pb(gammaResolution,sigmaResolution,nuResolution);
	pb.SetName(std::vector<std::string>{"  ","  ","  "});
	Log("\tInitialising the Probability Array\n");
	ErroredBinomial EB(d.maxK,res,bounder,qMax,qMax*nuMax,settings.ParallelWorkers);
	std::vector<double> nus;
	std::vector<std::vector<double>> searchVectors;
	Log("\tBeginning the search\n");
	for (int gamIdx = 0; gamIdx < gammaResolution; ++gamIdx)
	{
		double gamma = (gammaMin + deltaGamma * gamIdx);
		for (int sigIdx = 0; sigIdx < sigmaResolution; ++sigIdx)
		{
			double sigma = (sigmaMin + deltaSig * sigIdx);
			std::vector<double> probs;
			
			EB.Populate(gamma,sigma,dist); //generates a new probability array
			for (double nuIdx = 0; nuIdx < nuResolution; ++nuIdx)
			{
				double nu = nuMin + nuIdx * deltaNu;
				auto assign = AssignQs(nu,d,EB,alpha,L,accelerator,dist);
				// std::cout << nu << "  " << sigma << " scored " << assign.Score << std::endl;
				assign.Nu = nu;
				assign.Gamma = gamma;
				if (sigIdx == 0)
				{
					nus.push_back(nu);
				}
				probs.push_back(abs(assign.Score));

				if (assign.Score > best.Score)
				{
					best = assign;
				}
				if (globalVerbose)
				{
					pb.Update(gamIdx,sigIdx,nuIdx);
				}
			}
			searchVectors.push_back(probs);
		}
	}
	logPloidyPrior = 0;
	double bGamma = best.Gamma;
	Log("\tHarmonic search complete\n\tIdeal parameters are nu=" << best.Nu << ", sigma=" << best.Sigma << " gamma=" << best.Gamma << "\nBeginning high-resolution search\n");
	ErroredBinomial EB2(d.maxK,3*res,bounder,qMax,qMax*nuMax,settings.ParallelWorkers);
	EB2.Populate(best.Gamma,best.Sigma,dist);
	Log("\tHigh-resolution probability grid generated\n")
	best = AssignQs(best.Nu,d,EB2,alpha,L,1,dist);
	best.Gamma = bGamma;
	Log("\tCompleted high-resolution search")
	
	// for (int c = 0; c < d.Chromosomes.size(); ++c)
	// {
	// 	std::cout << "Break list for chrom " << d.Chromosomes[c].Name << std::endl;
	// 	chr_int prev = 0;
	// 	int prevQ = 0;
	// 	for (int i = 0; i < best.List[c].Index.size(); ++i)
	// 	{
	// 		chr_int p = best.List[c].Index[i];
	// 		std::cout << "\t" << prevQ << "->" << best.List[c].Resonance[i] << " at " << p << "(+" << p - prev << ")" << std::endl;
	// 		prev = p;
	// 		prevQ = best.List[c].Resonance[i];
	// 	}
	// }

	Log("Saving files to output\n");
	//check if output needs a dir creating
	auto dirSplit = JSL::split(settings.OutputName,'/');
	if (dirSplit.size() > 1)
	{
		std::string frontLoad = "";
		for (int i =0; i < dirSplit.size()-1; ++i)
		{
			frontLoad += dirSplit[i] + "/";
		}
		JSL::mkdir(frontLoad);
	}

	std::string metaFile = settings.OutputName + ".meta";
	JSL::initialiseFile(metaFile);
	std::ostringstream buffer;
	buffer << std::setprecision(8);
	for (int j = 0; j < searchVectors.size(); ++j)
	{
		for (int i =0; i < nus.size(); ++i)
		{
			buffer << nus[i] << " " << (sigmaMin + deltaSig * j);;
		
			buffer << " " << searchVectors[j][i];
			buffer << "\n";
		}
	}
	JSL::writeStringToFile(metaFile,buffer.str());

	std::string treeFile = settings.OutputName + ".tree";
	buffer.str("");
	buffer.clear();
	buffer << "Data " << settings.DataFile << "\n";
	std::vector<std::string> params = {"Nu","Sigma","Gamma","L"};
	std::vector<double> vals = {best.Nu, best.Sigma,best.Gamma,(double)L};
	for (int i = 0; i < params.size(); ++i)
	{
		buffer << params[i] << " " << vals[i] << "\n";
	}
	
	buffer << "===================\n";
	for (int c = 0; c < best.List.size(); ++c)
	{
		double prev = 0;
		for (int i = 0; i < best.List[c].Index.size(); ++i)
		{
			prev = best.List[c].Resonance[i];
			buffer << d.Chromosomes[c].Name << " " << best.List[c].Index[i] << " " << prev << "\n";
		}
		buffer << d.Chromosomes[c].Name << " " << d.Chromosomes[c].maxIdx << " " << prev << "\n";
	}
	JSL::initialiseFile(treeFile);
	JSL::writeStringToFile(treeFile,buffer.str());
	Log("Harmonic Analysis Complete\n")
}