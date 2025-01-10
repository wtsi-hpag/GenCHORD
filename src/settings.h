#pragma once
#include "JSL.h"
/*
	IMPORTANT! 
	Settings are defined programmatically in settings.dat, and this file is then constructed using the settingBuilder.py file. 
	This avoids boilerplate bloat, and integrates command line settings with the help files automatically. 
	DO NOT EDIT THIS FILE!
*/
#define SETTINGS_LIST \
	Setting(int,Ploidy,2,"ploidy")  /*The (assumed) default ploidy of the genome.*/\
	Setting(std::string,DataFile,"_no_file_","file")  /*The input datafile to be read*/\
	Setting(int,AccumulationFactor,100,"accumulate")  /*The size of the Accumulation Window*/\
	Setting(std::string,OutputDirectory,"Output","output")\
	Setting(bool,CreateArchive,1,"archive")  /*Controls if an archive file is created*/\

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