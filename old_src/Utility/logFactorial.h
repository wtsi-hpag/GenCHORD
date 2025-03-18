#pragma once
#include <vector>
#include <math.h>
class LogFactorial
{
	public:
		LogFactorial()
		{
			Populate(10);
		}
		double operator () (int k)
		{
			if (k > CurrentMax)
			{
				// std::cout << "Adding to buffer" << std::endl;
				Populate(k+10);
			}
			return Data[k];
		}
	private:
		int CurrentMax = 0;
		std::vector<double> Data;	
		void Populate(int end)
		{
			double prev = 0;
			if (CurrentMax > 0)
			{
				prev = Data[CurrentMax];
			}
			Data.resize(end+1,0);
			for (int i = CurrentMax+1; i <= end; ++i)
			{
				prev += log(i);
				Data[i] = prev;
			}
			CurrentMax = end;
		}
};