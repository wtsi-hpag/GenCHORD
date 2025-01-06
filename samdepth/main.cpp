#include <iostream>
#include "JSL.h"

std::vector<int> getHops(int argc, char ** argv)
{
	std::vector<int> hops = {10,100,1000};
	for (int k = 1; k < argc; ++k)
	{
		try
		{
			int q = std::stoi(argv[k]);
			if (JSL::FindXInY(q,hops) == -1)
			{
				hops.push_back(q);
			}
		}
		catch (...)
		{
			//don't
		}
	}
	return hops;
}

struct Hopper
{
	std::ofstream File;
	int Counter;
	int Sum;
	int HopSize;
	long long int Position;

	Hopper(int hop)
	{
		HopSize = hop;
	}
	void NewFile(std::string fileName)
	{
		File.close();
		File.open(fileName);
		Counter = 0;
		Sum = 0;
		Position = 0;
	}
	void Close()
	{
		File.close();
	}
	void Update(int idx, int k)
	{
		if (Position == 0)
		{
			Position = idx;
		}
		Sum += k;
		Counter += 1;
		if (Counter == HopSize)
		{
			File << Position << " " << Sum << "\n";
			Counter = 0;
			Sum = 0;
			Position = 0;
		}
	}
};

int main(int argc, char ** argv)
{
	
	if (!JSL::PipedInputFound())
	{
		std::cout << "Please pipe in the output from a samtools depth command" << std::endl;
		return -1;
	}

	JSL::Argument<std::string> Target("output","o",argc,argv);
	auto hops = getHops(argc,argv);

	JSL::mkdir(Target.Value);
	std::vector<Hopper> crawler;
	for (auto hop : hops)
	{
		crawler.push_back(Hopper(hop));
	}

	std::string prev = "";
	forLineInPipedInput(
		std::vector<std::string> vec = JSL::split(PIPE_LINE,' ');
		// if (vec.size() != 3)
		// {

		// }
		if (vec[0] != prev)
		{
			prev= vec[0];

			for (int i = 0; i < crawler.size(); ++i)
			{
				std::string name = Target.Value + "/" + prev + "_" + std::to_string(hops[i]) + ".dat";
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

	//the complicated archiving bit
	auto split = Target.Value.find_last_of('/');
	std::string directory = Target.Value.substr(0,split+1);
	std::string name = Target.Value.substr(split+1);
	std::string dirString = " ";
	if (directory.size() > 0)
	{
		dirString = " -C " + directory + " ";
	}
	std::string cmd = "tar -cf " + Target.Value + ".gcd" + dirString + name + "; rm -r " + Target.Value;

	system(cmd.c_str());
		
	
	return 0;
	
}