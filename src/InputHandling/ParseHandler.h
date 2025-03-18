#pragma once
#include <vector>
#include <sstream>
#include <regex>
#include "../settings.h"
#include "../Utility/StringSanitiser.h"
#include "JSL.h"
#include "Data.h"
#include "RawFileParser.h"

//this function detects piped input or assigned files, and then reads them, parses and packages them into the data entity.
DataHolder ParseData();