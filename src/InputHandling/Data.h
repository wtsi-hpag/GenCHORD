#pragma once

#include <tuple>
#include <vector>
#include "../settings.h"
#include "../Utility/Log.h"

typedef unsigned long long int lint;
typedef lint dnaindex; 

template<class Xtype, class Ytype> 
struct XY
{
	std::vector<Xtype> X;
	std::vector<Ytype> Y;
	XY(int size=0){X.resize(size); Y.resize(size);}
	void Add(Xtype x, Ytype y){X.push_back(x); Y.push_back(y);}
	void Resize(int n){X.resize(n); Y.resize(n);}
};
class Datum
{
	public:
		lint Coverage;
		lint SquareSum;
		dnaindex Index;
		Datum(): Coverage(0), Index(0), SquareSum(0){};
		Datum(dnaindex idx, lint k): Coverage(k), Index(idx), SquareSum(pow(k,2)){};
		Datum(dnaindex idx, lint k, lint kSq) : Index(idx), Coverage(k), SquareSum(kSq){};
		void Set(dnaindex idx, lint k)
		{
			Index = idx;
			Coverage = k;
			SquareSum = pow(k,2);
		}
};

class CoverageArray
{
	public:
		std::string Name;
		CoverageArray(){};
		CoverageArray(const std::string & name) : Name(name){Data.resize(0);};
		CoverageArray(const std::string & name, const std::vector<std::tuple<dnaindex, lint, lint>> & data);
		Datum operator [](int i) const {return Data[i];};
		Datum & operator[](int i) {return Data[i];};
		void AddData(dnaindex Index,unsigned int k);
		int size() const;
		void Statistics();
		double AggregateMean;
		double AggregateVariance;
		double RawMean;
		double RawVariance;
		void FlagTruncated();
		XY<lint,lint> GetCoverage(int skipFactor=1);
	private:
		std::vector<Datum> Data;

};



class DataHolder
{
	public: 
		DataHolder();
		DataHolder(const std::vector<CoverageArray> & data);
		void NewChromosome();
		void Append(const std::string & name, const std::vector<std::tuple<dnaindex,lint,lint>> & element);	
		CoverageArray operator[](int i) const {return data[i];};
		CoverageArray & operator[](int i){return data[i];};
		std::vector<int> Histogram() const;
		std::vector<int> Histogram(int chromosome) const;
		size_t size();
		void Analyse();
		void TruncateHistogram(std::vector<int> & vec, lint nTotal) const;
		double OverallMean;
	private:
		std::vector<CoverageArray> data;
		void internalHistogram(int c, std::vector<int> & output) const;
		bool TruncateMode = true;
};

// typedef std::vector<CoverageArray> DataHolder;
