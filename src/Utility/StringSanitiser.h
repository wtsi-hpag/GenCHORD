#pragma once
#include <iostream>
#include <string>
#include <regex>

inline bool StringIsSanitised(const std::string & inputString)
{
	 std::regex unsafePattern("[;&|<>$`\\\"\\'\\\\]");
    
    // Use std::regex_search to find any unsafe characters
    return !std::regex_search(inputString, unsafePattern);
}