#pragma once
#include <vector>
#include <string>
#include "../Utility/basicFunctions.h"
#include "../../libs/JSL/JSL.h"

#include "../settings.h"
#include "ChromosomeCoverage.h"

struct FileReading
{
	long long int CurrentLine;
	int CurrentChromosome;
	std::string ChromSignifier; //chromosomes not always labelled using integers

	int GapSize;
	chr_int PreviousIndex;
	bool ChromActive;

	double CoverageSum;
	double CoverageSquareSum;
	double SmoothingAverage;
	long long int LoadedData;
	long long int SpoofedData;
	int MaxCoverage;
	// long long int
	ChromosomeCoverage Chromosome;
};

class Data
{
	public:
		double mean;
		std::vector<ChromosomeCoverage> Chromosomes;
		int maxK=0;
		double Mean;
		double Deviation;
		bool Loaded;
		Data();

		Data(Settings & settings);

		Data(std::string target, int thinning, std::string targetChromosome,double smoothing);

	private:

		void PrepareForLoading();

		void ParseLine(const std::vector<std::string> & line, const Settings & settings);

		FileReading File;

		void AddData(int index, int coverage, const Settings & settings, bool spoofed);

};	