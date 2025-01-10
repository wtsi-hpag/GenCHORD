#include <iostream>
#include "JSL.h"
#include "Utility/Log.h"
#include "settings.h"


void ConfigureLogging()
{
	LOGCFG.headers = Settings.LogHeaders;
	LOGCFG.SetLevel(Settings.LogLevel);
}

int main(int argc, char ** argv)
{
	Settings.Configure(argc,argv);
	ConfigureLogging();

	
	std::cout << Settings.Ploidy << " " << Settings.OutputDirectory << " " << std::endl;
	return 0;
}