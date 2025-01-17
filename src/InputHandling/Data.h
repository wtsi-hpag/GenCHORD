#pragma once
#include <vector>
#include <vector>

#include "../settings.h"
#include "../Utility/Log.h"

typedef unsigned long long int lint;
typedef unsigned long long int dnaindex; 

class Datum
{
	public:
		lint Coverage;
		lint SquareSum;
		dnaindex Index;
		Datum(): Coverage(0), Index(0), SquareSum(0){};
		Datum(dnaindex idx, lint k): Coverage(k), Index(idx), SquareSum(pow(k,2)){};
		Datum(dnaindex idx, lint k, lint kSq) : Index(idx), Coverage(k), SquareSum(kSq){};
};

class CoverageArray
{
	public:
		std::string Name;
		CoverageArray(){};
		CoverageArray(const std::string & name) : Name(name){Data.resize(0);};
		CoverageArray(const std::string & name, std::vector<std::tuple<dnaindex, lint, lint>> & data);
		Datum operator [](int i) const {return Data[i];};
		Datum & operator[](int i) {return Data[i];};
		void AddData(dnaindex Index,unsigned int k);
		int size();
		void Statistics();
		double AggregateMean;
		double AggregateVariance;
		double RawMean;
		double RawVariance;
		void FlagTruncated();
	private:
		std::vector<Datum> Data;

};



class DataHolder
{
	public: 
		DataHolder();
		DataHolder(const std::vector<CoverageArray> & data);
		void Append(const std::string & name, std::vector<std::tuple<dnaindex,lint,lint>> & element);	
		CoverageArray operator[](int i) const {return data[i];};
		CoverageArray & operator[](int i){return data[i];};
		std::vector<int> Histogram();
		std::vector<int> Histogram(int chromosome);

		void Analyse();
	private:
		std::vector<CoverageArray> data;
		void internalHistogram(int c, std::vector<int> & output);
};

// typedef std::vector<CoverageArray> DataHolder;
