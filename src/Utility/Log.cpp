#include "Log.h"

void structlog::SetLevel(int i)
{
	LOG(INFO) << "Testing " << i << "\n";
	switch(i){
		case 0: 
			level=ERROR; break;
		case 1:
			level=WARN;break;
		case 2:
			level=INFO;break;
		case 3:
			level=DEBUG;break;
		default:
			LOG(WARN) << "\"" << i << "\" is not a valid logging level (0-3). Default to INFO.\n";
			level=INFO;break;
	}
}

structlog LOGCFG = {};

