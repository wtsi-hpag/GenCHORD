#include "Data.h"

void CoverageArray::AddData(dnaindex idx, unsigned int k)
{
	Data.push_back(Datum(idx,k));
}
void CoverageArray::FlagTruncated()
{
	// Data.pop_back();
}

CoverageArray::CoverageArray(const std::string & name, const std::vector<std::tuple<dnaindex, int, int>> & data) : Name(name)
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
double CoverageArray::GetMean()
{
	double mu = 0;
	double muSq = 0;
	for (int i = 0; i < Data.size(); ++i)
	{
		mu += Data[i].Coverage;
		muSq += Data[i].SquareSum;
	}
	// Mean = mu/Data.size();
	// double trueMean = Mean/Settings.AccumulationFactor;
	// Variance = muSq/(Data.size()*Settings.AccumulationFactor) - trueMean * trueMean;
	LOG(DEBUG) << "Sum " << mu << " sqSum " << muSq << " var " << muSq/(Data.size()) - pow(mu/Data.size(),2);
	return Mean;

}

DataHolder::DataHolder()
{
	data.resize(0);
}
DataHolder::DataHolder(const std::vector<CoverageArray> & input) : data(input)
{

}
void DataHolder::Append(const std::string & name, const std::vector<std::tuple<dnaindex,int,int>> & element)
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
	double muTotal = 0;
	double sqTotal = 0;
	double longTotal = 0;
	for (int c = 0; c < data.size(); ++c)
	{
		data[c].GetMean();
		// LOG(INFO) << "Chromosome " << data[c].Name << " has aggregate " << data[c].Mean << "±" << sqrt(data[c].Variance) << " corresponding to mean coverage " << data[c].Mean / Settings.AccumulationFactor;
		// int n = data[c].size();
		// muTotal += data[c].Mean * n;
		// longTotal += n;
		// sqTotal += (data[c].Variance + data[c].Mean*data[c].Mean)*n;
	}


}