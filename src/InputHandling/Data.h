#pragma once
#include <vector>
typedef unsigned long long int dnaindex; 
class Datum
{
	public:
		unsigned int Coverage;
		unsigned int SquareSum;
		dnaindex Index;
		Datum(): Coverage(0), Index(0), SquareSum(0){};
		Datum(int idx, int k): Coverage(k), Index(idx), SquareSum(k*k){}
};

class CoverageArray
{
	public:
		Datum operator [](int i) const {return Data[i];};
		Datum & operator[](int i) {return Data[i];};
		void AddData(dnaindex Index,unsigned int k);
		// void AccumulateData(dna)
	private:
		std::vector<Datum> Data;

};