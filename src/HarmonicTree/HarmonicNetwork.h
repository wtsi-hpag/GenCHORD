#pragma once
#include <vector>
#include "Path.h"
#include "../data.h"

class HarmonicNetwork
{
	public:
		HarmonicNetwork(){};
		HarmonicNetwork(const Data & d, int assignedChromosome, int qMax, int L, bool scan);

		void Initialise(const Data & d, int assignedChromosome, int qMax, int L,bool scan);



		void Navigate(const Data & d, double nu, double gamma);

		Path BestPath();
		bool ScanMode;
	private:
		std::vector<std::vector<Path>> Paths;

		int Memory;
		int Qmax;
		int BufferSize;
		chr_int DataSize;
		int MyChromosome;

};