#include <iostream>
#include <stdexcept>
#include "JSL.h"
#include "Aggregator.h"
#include <filesystem>
std::vector<int> getWindows(int argc, char ** argv)
{
	std::vector<int> windows = {100,1000,10000};
	for (int k = 1; k < argc; ++k)
	{
		try
		{
			int q = std::stoi(argv[k]);
			if (JSL::FindXInY(q,windows) == -1)
			{
				windows.push_back(q);
			}
		}
		catch (...)
		{
			//don't do anything -- errors will most likely be arising from other command line arguments 
		}
	}
	return windows;
}


int main(int argc, char ** argv)
{
	if (!JSL::PipedInputFound())
	{
		std::cerr << "DEPTHSPLITTER ERROR: Please pipe in the output from a samtools depth command" << std::endl;
		return 1;
	}
	try
	{
		JSL::Argument<std::string> Target("output","o",argc,argv);
		JSL::Argument<std::string> Delimiter("\t","delim",argc,argv);

		auto windows = getWindows(argc,argv);

		//the complicated archiving bit
		auto split = Target.Value.find_last_of('/');
		std::string directory = Target.Value.substr(0,split+1);
		std::string name = Target.Value.substr(split+1);
		if (directory.size() > 0)
		{
			JSL::mkdir(directory);
		}

		std::ofstream tar(Target.Value + ".gcd",std::ios::binary);
		std::vector<Aggregator> crawler;

		for (auto window : windows)
		{
			crawler.push_back(Aggregator(window,tar));
		}
		std::string prev = "";

		char delim = Delimiter.Value[0];
		forLineInPipedInput(
			std::vector<std::string> vec = JSL::split(PIPE_LINE,delim);
			if (vec.size() != 3)
			{
				throw std::runtime_error("Expected samtools depth format of tab-delimited 3-column output. You likely got the wrong delimiter");
			}
			if (vec[0] != prev)
			{
				prev= vec[0];

				for (int i = 0; i < crawler.size(); ++i)
				{
					std::string name = prev + "_" + std::to_string(windows[i]) + ".dat";
					crawler[i].NewFile(name);
				}
			}
			int idx = std::stoi(vec[1]);
			int k = std::stoi(vec[2]);
			for (int i = 0; i < crawler.size(); ++i)
			{
				crawler[i].Update(idx,k);
			}

		);
	

	

		for (int i = 0; i < crawler.size(); ++i)
		{
			crawler[i].Close();
		}	

		tar_to_stream_tail(tar);
		
		return 0;
	}
	catch (const std::exception& e) 
	{
		std::cerr << "ERROR: " << e.what() << std::endl;
		return 1;
	}
}
