#include <iostream>
#include "JSL.h"
#include "Utility/Log.h"
#include "settings.h"


int main(int argc, char ** argv)
{
	// process settings

	LOGCFG.headers = true; 
	LOGCFG.level = DEBUG;

	// Settings settings(argc,argv);
	Settings.Configure(argc,argv);

	std::cout << Settings.Ploidy << " " << Settings.OutputDirectory << " " << std::endl;
	return 0;
}