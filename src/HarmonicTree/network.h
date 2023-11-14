#pragma once
#include "node.h"
#include "../data.h"


class Network
{
	public:
		Network(){};
		Network(const Data & d, int chrom, int qMax)
		{
			Init(d,chrom,qMax);
		}

		void Init(const Data & d, int chrom, int qMax)
		{
			Nlayers = d.Chromosomes[chrom].Counts.size();
			Nodes.resize(Nlayers);
			for (int i = 0; i < Nlayers; ++i)
			{
				Nodes[i].resize(qMax);
				for (int q = 0; q < qMax; ++q)
				{
					Nodes[i][q].Assign(q,i);
				}
			}
		}
	private:
		int Nlayers;
		std::vector<std::vector<Node>> Nodes;
		
};