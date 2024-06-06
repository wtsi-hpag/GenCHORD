#pragma once
#include <vector>
#include <iostream>
#include <sstream>
struct Coordinate
{
	int Index;
	short Value;
	Coordinate(){};
	Coordinate(int idx, short val) : Index(idx), Value(val)
	{
		
	}
};

class Path
{
	public:
		double Score;
		double CumulativeLinearScore;
		int CurrentQ;
		int Idx;
		int CoordSum = -1;
		int QSum = -1;
		double Nu;
		std::vector<Coordinate> Route;
		// std::vector<int> PreviousIdx;
		Path()
		{
			Score = 0;
			CumulativeLinearScore = 0;
			CurrentQ = -1;
		}

		void InitialStep(double score, int q)
		{
			Score = score;
			CumulativeLinearScore = score;
			CurrentQ = q;
			Route = {Coordinate(0,q)};
			// PreviousIdx.resize(0);
			// PushJump(q,0);
			Idx = 0;
			CoordSum = 0;
			QSum = q;
		}

		void RecordScore(double score, int q, int idx)
		{
			Score = score;
			Idx = idx;
			CurrentQ =q;
		}

		void AddStep(double score, int q, int idx,const Path * node,int specIdx)
		{
			Score = score;
			CurrentQ = q;
			Idx = idx;
			bool jumped = q != node->CurrentQ;
			if (jumped || (CoordSum != node->CoordSum) || (QSum != node->QSum))
			{
				Route = node->Route;
				CoordSum = node->CoordSum;
				QSum = node->QSum;
			}

			if (jumped)
			{
				CoordSum += specIdx;
				QSum += q;
				Route.push_back(Coordinate(specIdx,q));
			}
		}

		void PushJump(int q, int idx)
		{
			
			
			// PreviousIdx.push_back(idx);
		}

		std::string TreeOutput(std::string chrName)
		{
			std::ostringstream os;

			for (int i = 0; i < Route.size(); ++i)
			{
				os << chrName << " " << Route[i].Index << " " << Route[i].Value << "\n";
			}
			return os.str();
	

		}
};