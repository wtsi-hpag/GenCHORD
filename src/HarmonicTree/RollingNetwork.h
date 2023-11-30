#pragma once
#include <vector>
#include "Path.h"
#include "../data.h"

class RollingNetwork
{
	public:
		RollingNetwork(){};
		RollingNetwork(const Data & d, int assignedChromosome, int qMax, int L);

		void Initialise(const Data & d, int assignedChromosome, int qMax, int L);



		void Navigate(const Data & d, double nu, double gamma);

		Path BestPath();
	private:
		std::vector<std::vector<Path>> Paths;

		int Memory;
		int Qmax;
		int BufferSize;
		chr_int DataSize;
		int MyChromosome;

};