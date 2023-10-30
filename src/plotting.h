#pragma once
#include "JSL.h"
#include "data.h"
#include "Transitions.h"
#include "settings.h"
extern int nPlot;
void basicPlot(JSL::gnuplot & gp,Data & d,int chrom);

void TransitionPlot(JSL::gnuplot & gp,std::vector<chr_int> idx, std::vector<double> freq,double nuCorrect,std::string legend);

void recoverMeta(Settings & settings);

class TreeMeta
{
	public:
		std::string DataFile;
		double Nu;
		
		double Sigma;
		double Gamma;
		int L;

		TreeMeta(Settings & s);
		std::vector<std::vector<int>> Harmonic;
		std::vector<std::vector<double>> Frequency;
		std::vector<std::vector<chr_int>> Idx;
		std::vector<std::string> Name;
		std::vector<int> TrueChrom;
	
};

void TransitionFrame(Data & d, Settings & settings, std::vector<chr_int> idx, std::vector<double> freq, int c,bool instantPlot);

void OutputPlot(Data & d,Settings & settings);

void ComparisonPlots(Settings & settings);