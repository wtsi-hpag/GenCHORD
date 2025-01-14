#pragma once
#include <vector>

#include "../Utility/Log.h"

typedef unsigned long long int dnaindex; 
class Datum
{
	public:
		unsigned int Coverage;
		unsigned int SquareSum;
		dnaindex Index;
		Datum(): Coverage(0), Index(0), SquareSum(0){};
		Datum(dnaindex idx, int k): Coverage(k), Index(idx), SquareSum(k*k){};
		Datum(dnaindex idx, int k, int kSq) : Index(idx), Coverage(k), SquareSum(kSq){};
};

class CoverageArray
{
	public:
		CoverageArray(){Data.resize(0);};
		CoverageArray(const std::vector<std::tuple<dnaindex, int, int>> & data);
		Datum operator [](int i) const {return Data[i];};
		Datum & operator[](int i) {return Data[i];};
		void AddData(dnaindex Index,unsigned int k);
		int Size();
		// void AccumulateData(dna)
	private:
		std::vector<Datum> Data;

};

typedef std::vector<CoverageArray> DataHolder;