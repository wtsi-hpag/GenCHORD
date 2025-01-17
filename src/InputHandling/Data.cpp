#include "Data.h"

void CoverageArray::AddData(dnaindex idx, unsigned int k)
{
	Data.push_back(Datum(idx,k));
}
void CoverageArray::FlagTruncated()
{
	Data.pop_back();
}

XY<lint, lint> CoverageArray::GetCoverage(int skipFactor)
{
	XY<lint, lint> out;
	// std::vector<lint> out(Data.size()/skipFactor);
	for (int i = 0; i < Data.size(); i+=skipFactor)
	{
		out.Add(Data[i].Index,Data[i].Coverage);
	}
	return out;
}

CoverageArray::CoverageArray(const std::string & name, const std::vector<std::tuple<dnaindex, lint, lint>> & data) : Name(name)
{
	Data.resize(data.size());
	for (int i = 0; i < data.size(); ++i)
	{
		Data[i] = Datum(std::get<0>(data[i]),std::get<1>(data[i]),std::get<2>(data[i]));
	}
}
int CoverageArray::size()
{
	return Data.size();
}
void CoverageArray::Statistics()
{
	double mu = 0;
	double aggSq = 0;
	double muSq = 0;
	int logVal = 2;
	double lv = log(logVal);
	std::vector<int> hist(100,0);
	lint maxK = 0;
	for (int i = 0; i < Data.size(); ++i)
	{
		lint k = Data[i].Coverage;
		if (k > maxK)
		{
			maxK = k;
		}

		// if (k >= hist.size())
		// {
		// 	LOG(DEBUG) << "increasing size to " << k;
		// 	hist.resize(k,0);
		// }
		// hist[k]+=1;
		// LOG(DEBUG) << hist[k];
		mu += k;
		aggSq += pow((double)k,2);
		muSq += Data[i].SquareSum;
		
	}
	AggregateMean = mu/Data.size();
	RawMean = AggregateMean / Settings.AccumulationFactor;

	AggregateVariance = aggSq/Data.size() - pow(AggregateMean,2);
	RawVariance = muSq/(Data.size() * Settings.AccumulationFactor) - pow(RawMean,2);
	
	// LOG(DEBUG) << "Beginning hist analysis " << maxK;
	// double target = Settings.TruncationFactor;
	// int cutoff = 0;
	// double cumsum = hist[cutoff] * 1.0/Data.size();
	// while (cutoff < hist.size() && cumsum < target)
	// {

	// 	cumsum += hist[cutoff] * 1.0/Data.size();
	// 	cutoff +=1;
	// }
	// hist.resize(cutoff);
	
	LOG(INFO) << "Chr: "<<Name << "\tAggregate: " << std::setprecision(5) << AggregateMean<< " ± " << std::setprecision(4) << sqrt(AggregateVariance)/AggregateMean * 100 << "%\tRaw: " << std::setprecision(5) << RawMean << " ± " << std::setprecision(4) << sqrt(RawVariance)/RawMean*100 << "%";
}

DataHolder::DataHolder()
{
	data.resize(0);
}
DataHolder::DataHolder(const std::vector<CoverageArray> & input) : data(input)
{

}
void DataHolder::Append(const std::string & name, const std::vector<std::tuple<dnaindex,lint,lint>> & element)
{
	data.emplace_back(CoverageArray(name, element));
}

void DataHolder::internalHistogram(int c, std::vector<int> & output)
{
	for (int i = 0; i < data[c].size(); ++i)
	{
		int k = data[c][i].Coverage;
		if (k >= output.size())
		{
			output.resize(k+1,0);
		}
		output[k] +=1;
	}
}
std::vector<int> DataHolder::Histogram(int c)
{
	std::vector<int> output;
	internalHistogram(c,output);
	return output;
}
std::vector<int> DataHolder::Histogram()
{
	std::vector<int> output;
	for (int c =0 ; c< data.size(); ++c)
	{
		internalHistogram(c,output);
	}
	return output;
}

void DataHolder::Analyse()
{
	int c = 0;
	while (c < data.size())
	{
		if (data[c].size() > 1)
		{
			++c;
		}
		else
		{
			LOG(WARN) << "Chromosome " << data[c].Name << " has insufficient associated with it and has been removed.\n\tThis is not an error if you expect it to be shorter than twice the accumulation factor ("<< Settings.AccumulationFactor << "). \n\tOtherwise this is an indicator of a malformed datafile.";
			data.erase(data.begin()+c);
		}
	}
	double muTotal = 0;
	double sqTotal = 0;
	double longTotal = 0;
	for (int c = 0; c < data.size(); ++c)
	{
		LOG(DEBUG) << "Running stats for" << data[c].Name;
		data[c].Statistics();
		
		// int n = data[c].size();
		// muTotal += data[c].Mean * n;
		// longTotal += n;
		// sqTotal += (data[c].Variance + data[c].Mean*data[c].Mean)*n;
	}


}

size_t DataHolder::size()
{
	return data.size();
}