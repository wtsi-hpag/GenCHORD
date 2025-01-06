#include <iostream>
#include "JSL.h"
int main(int argc, char ** argv)
{
	JSL::Argument<std::string> Input("Data/Aaron.dat","-f");
	std::string command = "cat " + Input.Value + " | ./depthsplitter";

	system(command.c_str());

	return 0;
}