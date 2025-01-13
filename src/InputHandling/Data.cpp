#include "Data.h"

void CoverageArray::AddData(dnaindex idx, unsigned int k)
{
	Data.push_back(Datum(idx,k));
}

CoverageArray::CoverageArray(const std::vector<std::tuple<dnaindex, int, int>> & data)
{
	Data.resize(data.size());
	for (int i = 0; i < data.size(); ++i)
	{
		Data[i] = Datum(std::get<0>(data[i]),std::get<1>(data[i]),std::get<2>(data[i]));
	}
}