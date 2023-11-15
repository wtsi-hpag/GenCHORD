#pragma once

class Node
{
	public:
		double Score;
		Node * Prev;
		int Idx;
		double NodeProb;
		double CumulativeNodeProb=0;
		short int Q;
		// bool Connected = false;
		Node(){}
		
		void Assign(int q, int id)
		{
			Q = q;
			Idx = id;
			Score = 0;
		}

		
		// void UpdateProb()
};