#include "Data.h"

void CoverageArray::AddData(dnaindex idx, unsigned int k)
{
	Data.push_back(Datum(idx,k));
}

CoverageArray::CoverageArray(const std::vector<std::tuple<dnaindex, int, int>> & data)
{
	// LOG(DEBUG) << "I Should be size" << data.size();
	Data.resize(data.size());
	for (int i = 0; i < data.size(); ++i)
	{
		Data[i] = Datum(std::get<0>(data[i]),std::get<1>(data[i]),std::get<2>(data[i]));
	}
}
int CoverageArray::Size()
{
	return Data.size();
}