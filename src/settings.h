#pragma once
#include "JSL.h"
#include "Utility/Log.h"
/*
	IMPORTANT! 
	
	Settings are defined using some funky macros which remove the need to declare/initialise variables on multiple lines.

	To add a new setting variable, add it in following the pattern Setting(type,VariableName,default value,command-line-string. Everything will be sorted such that a member variable of Settings will be constructed with name Variable name, of type "type" and initialised to the default value. 

	If you pass the command-line-string as a command line argument (or in a config file), the Settings object will instead use that value, using the most recent value it encounters (left-right in CLA, top-down in config file). 
*/
#define SETTINGS_LIST \
	Setting(int,Ploidy,2,"ploidy")  /*The (assumed) default ploidy of the genome.*/\
	Setting(std::string,DataFile,"_no_file_","file")  /*The input datafile to be read*/\
	Setting(int,AccumulationFactor,100,"accumulate")  /*The size of the Accumulation Window*/\
	Setting(std::string,Output,"genchord_output","output")\
	Setting(bool,CreateArchive,1,"archive")  /*Controls if an archive file is created*/\
	Setting(int,LogLevel,2,"log")  /*The logging level (0=errors only, 1 = errors + warnings, 2=general info, 3=debug mode)*/\
	Setting(char,StreamDelimiter,'\t',"delim")	\
	Setting(bool,LogHeaders,1,"loghead")  /*Controls if the logging headers are printed*/\
	Setting(int,PlotBinFactor,1,"plotbin") \
	Setting(double,TruncationFactor,0.999,"truncate")/*Controls the amount of high-coverage tail which is truncated*/\
	Setting(int,DataGap,1,"datagap")/*the expected gap between genome indices in the input file*/\
	Setting(std::string,IgnoreChromosomeFlag,"CACP","ignore")\
	Setting(int,AutocorrelationLength,10,"autocorr")/*The (approximate) lengthscale over which the coverage data are not independent*/\
	Setting(double, ErrorMax,1e-3,"epsilon") /*the maximum permitted error rate during inference*/ \
	Setting(double, ContaminationMin,-0.8,"cont-min") /* the lower deviation of contamination*/ \
	Setting(double, ContaminationMax,0.8,"cont-max") /*The upper deviation of contamination*/\
	Setting(bool,ProcessMode,false,"process") /*activates a special mode which only reads in files*/\
	Setting(bool,TreeMode,true,"tree") \
// Do not edit anything below this line!
class SettingsObject		
{
	public:
		#define Setting(type, name, defaultValue, cmdArg) type name;
			SETTINGS_LIST
		#undef Setting

		// Constructor
		SettingsObject(){};
		SettingsObject(int argc, char** argv)
		{
			Configure(argc,argv);
		}
		
		void Configure(int argc, char**argv)
		{
			std::string configFile = JSL::Argument<std::string>("__none__","config",argc,argv);
			if (configFile.size() == 0 || configFile == "__none__")
			{
				Initialise(argc, argv);
			}
			else
			{
				Initialise(configFile,' ');
			}
		}

		template<class T, class U>
		void Initialise(T a, U b)
		{
			#define Setting(type, name, defaultValue, cmdArg) name = JSL::Argument<type>(defaultValue, cmdArg, a,b);
			SETTINGS_LIST
			#undef Setting
		}
};

extern SettingsObject Settings;