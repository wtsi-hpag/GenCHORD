#include "Log.h"

void structlog::SetLevel(int i)
{
	LOG(INFO) << "Testing " << i << "\n";
	switch(i){
		case 0: 
			level=ERROR;
		case 1:
			level=WARN;
		case 2:
			level=INFO;
		case 3:
			level=DEBUG;
		default:
			LOG(WARN) << "\"" << i << "\" is not a valid logging level (0-3). Default to INFO.\n";
			level=INFO;
	}
}

structlog LOGCFG = {};

