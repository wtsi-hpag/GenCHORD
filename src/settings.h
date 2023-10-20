#pragma once
#include "JSL.h"

extern bool globalRegress;
class Settings
{
	public:
		int Ploidy;
		double PloidyPrior;
		double sigmaMin;
		double sigmaMax;
		int sigmaResolution;
		int nuResolution;
		int Accelerator;
		double Gamma;
		int Qmax;
		int L;
		double alpha;

		std::string OutputName;
		std::string TargetChromosome;
		std::string DataFile;
		int DataThinning;
		bool Quiet;
		bool PlotOnly;
		std::string ComparePlot;
		double MemorySmoothing;
		int ParallelWorkers;
		Settings(){};

		void Initialise(int argc, char**argv)
		{
			Ploidy = JSL::Argument<int>(2,"ploidy",argc,argv);
			PloidyPrior= JSL::Argument<double>(0.5,"ploidyPrior",argc,argv);
			sigmaMin =JSL::Argument<double>(1,"sigmaMin",argc,argv);
			sigmaMax = JSL::Argument<double>(10,"sigmaMax",argc,argv);
			sigmaResolution= JSL::Argument<int>(3,"sigmaResolution",argc,argv);
			nuResolution = JSL::Argument<int>(100,"nuResolution",argc,argv);
			Accelerator = JSL::Argument<int>(1,"accelerate",argc,argv);
			Gamma =  JSL::Argument<double>(5,"gamma",argc,argv);
			Qmax = JSL::Argument<int>(10,"Qmax",argc,argv);
			L = JSL::Argument<int>(100000,"L",argc,argv);
			alpha = JSL::Argument<double>(1e-2,"alpha",argc,argv);
			TargetChromosome = JSL::Argument<std::string>("all","c",argc,argv);
			DataFile = JSL::Argument<std::string>("_no_file_provided_","f",argc,argv);
			DataThinning = JSL::Argument<int>(1,"thin",argc,argv);
			Quiet = JSL::Toggle(false,"q",argc,argv);
			OutputName = JSL::Argument<std::string>("output","o",argc,argv);
			ParallelWorkers = JSL::Argument<int>(0,"worker",argc,argv);
			PlotOnly = JSL::Toggle(false,"plot",argc,argv);
			ComparePlot = JSL::Argument<std::string>("__none__","compare",argc,argv);

			MemorySmoothing = JSL::Argument<double>(0,"smooth",argc,argv);
		}
		void Initialise(std::string configFile)
		{
			Ploidy = JSL::Argument<int>(2,"ploidy",configFile,' ');
			PloidyPrior= JSL::Argument<double>(0.5,"ploidyPrior",configFile,' ');
			sigmaMin =JSL::Argument<double>(1,"sigmaMin",configFile,' ');
			sigmaMax = JSL::Argument<double>(10,"sigmaMax",configFile,' ');
			sigmaResolution= JSL::Argument<int>(8,"sigmaResolution",configFile,' ');
			nuResolution = JSL::Argument<int>(100,"nuResolution",configFile,' ');
			Accelerator = JSL::Argument<int>(2,"accelerate",configFile,' ');
			Gamma =  JSL::Argument<double>(5,"gamma",configFile,' ');
			Qmax = JSL::Argument<int>(10,"Qmax",configFile,' ');
			L = JSL::Argument<int>(100000,"L",configFile,' ');
			alpha = JSL::Argument<double>(1e-2,"alpha",configFile,' ');
			TargetChromosome = JSL::Argument<std::string>("all","c",configFile,' ');
			DataFile = JSL::Argument<std::string>("_no_file_provided_","f",configFile,' ');
			DataThinning = JSL::Argument<int>(1,"thin",configFile,' ');
			Quiet = JSL::Argument<bool>(false,"q",configFile,' ');
			OutputName = JSL::Argument<std::string>("output","o",configFile, ' ');
			ParallelWorkers = JSL::Argument<int>(0,"worker",configFile,' ');
		}
};	