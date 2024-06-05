#pragma once
#include "../libs/JSL/JSL.h"

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

		double gammaMin;
		double gammaMax;
		int gammaResolution;
		int Accelerator;
		double Gamma;
		int Qmax;
		int L;
		double ContinuityPrior;

		std::string OutputName;
		std::string TargetChromosome;
		bool AllChromosomes;
		std::string DataFile;
		int DataThinning;
		bool Quiet;
		bool PlotOnly;
		std::string ComparePlot;
		double MemorySmoothing;
		int ParallelWorkers;
		char DataFileDelimiter;
		Settings(std::string configFile, int argc, char** argv)
		{
			if (configFile == "__none__")
			{
				Initialise(argc,argv);
			}
			else
			{
				Initialise(configFile);
			}
		}

		void Initialise(int argc, char**argv)
		{
			InitialiseParams(argc,argv);

			Quiet = JSL::Toggle(false,"q",argc,argv); //flags can only be used with the command-arg interface, not with config files
		}
		void Initialise(std::string configFile)
		{
			InitialiseParams(configFile,' ');
		}

	private:

			template<class T, class U>
			void InitialiseParams(T a, U b)
			{
				Ploidy = JSL::Argument<int>(2,"ploidy",a,b);
				PloidyPrior= JSL::Argument<double>(1,"ploidyPrior",a,b);
				sigmaMin =JSL::Argument<double>(1,"sigmaMin",a,b);
				sigmaMax = JSL::Argument<double>(10,"sigmaMax",a,b);
				sigmaResolution= JSL::Argument<int>(3,"sigmaResolution",a,b);
				nuResolution = JSL::Argument<int>(100,"nuResolution",a,b);
				gammaMin =JSL::Argument<double>(0.1,"gammaMin",a,b);
				gammaMax = JSL::Argument<double>(10,"gammaMax",a,b);
				gammaResolution= JSL::Argument<int>(5,"gammaResolution",a,b);
				Accelerator = JSL::Argument<int>(1,"accelerate",a,b);
				Gamma =  JSL::Argument<double>(-1,"gamma",a,b);
				Qmax = JSL::Argument<int>(10,"Qmax",a,b);
				L = JSL::Argument<int>(10000,"L",a,b);
				ContinuityPrior = JSL::Argument<double>(1e-2,"alpha",a,b);
				TargetChromosome = JSL::Argument<std::string>("all","c",a,b);

				AllChromosomes = (TargetChromosome == "all"); 

				DataFile = JSL::Argument<std::string>("_no_file_provided_","f",a,b);
				DataThinning = JSL::Argument<int>(1,"thin",a,b);
				
				OutputName = JSL::Argument<std::string>("output","o",a,b);
				ParallelWorkers = JSL::Argument<int>(0,"worker",a,b);
				ComparePlot = JSL::Argument<std::string>("__none__","compare",a,b);

				MemorySmoothing = JSL::Argument<double>(0,"smooth",a,b);

				Quiet = JSL::Argument<bool>(false,"quiet",a,b);

				DataFileDelimiter = JSL::Argument<char>(' ',"data-delimiter",a,b);
			}

			
			
};	