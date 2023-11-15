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
			QMax = qMax;
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

		Node *  GetBestPath()
		{
			int finalLayer = Nodes.size()-10;

			int bestScore = Nodes[finalLayer][0].Score;
			Node * bestNode = &Nodes[finalLayer][0];

			for (int q = 0; q < QMax; ++q)
			{
				std::cout << q << " has " << Nodes[finalLayer][q].Score << std::endl;
				if (Nodes[finalLayer][q].Score > bestScore)
				{
					bestScore = Nodes[finalLayer][q].Score;
					bestNode = & Nodes[finalLayer][q];
				}
			}
			return bestNode;
		}
	// private:
		int QMax;
		int Nlayers;
		std::vector<std::vector<Node>> Nodes;
		
};