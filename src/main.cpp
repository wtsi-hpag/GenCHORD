#include <iostream>
#include "JSL.h"
#include "Utility/Log.h"
#include "settings.h"


int main(int argc, char ** argv)
{
	// process settings

	LOGCFG.headers = true; 
	LOGCFG.level = DEBUG;

	Settings settings(argc,argv);

	std::cout << settings.Ploidy << " " << settings.OutputDirectory << " " << settings.Test << std::endl;
	return 0;
}