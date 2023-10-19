#pragma once
#include "JSL.h"
#include "data.h"
#include "Transitions.h"
#include "settings.h"

extern int nPlot;
inline void basicPlot(JSL::gnuplot & gp,Data & d,int chrom)
{
	gp.Plot(d.Chromosomes[chrom].Idx,d.Chromosomes[chrom].Counts,JSL::LineProperties::Legend("Data"),JSL::LineProperties::Colour("black"));
	
}

inline void TransitionPlot(JSL::gnuplot & gp,std::vector<chr_int> idx, std::vector<double> freq,int maxIdx)
{
	std::vector<chr_int> remap_idx;
	std::vector<double> remap_freq;
	double prev = freq[0];
	for (int i = 0; i < idx.size(); ++i)
	{
		remap_idx.push_back(idx[i]);
		remap_idx.push_back(idx[i]);
		remap_freq.push_back(prev);
		prev = freq[i];
		remap_freq.push_back(prev);
	}
	remap_idx.push_back(maxIdx);
	remap_freq.push_back(prev);
	gp.Plot(remap_idx,remap_freq,JSL::LineProperties::Legend("Treefit"));
}

inline void recoverMeta(Settings & settings)
{
	auto q = JSL::split(settings.DataFile,'.');
	std::string reconstructFile = q[0] + ".tree";
	settings.OutputName = q[0];
	int i = 0;
	forLineVectorIn(reconstructFile, ' ',

		if (FILE_LINE_VECTOR[0] == "Data")
		{
			settings.DataFile = FILE_LINE_VECTOR[1];
			return;
		}

		std::cout << "I could not find the metadata for the raw data, and so I am unable to plot" << std::endl;
		exit(1);
	);
}

inline void TransitionFrame(Data & d, Settings & settings, std::vector<chr_int> idx, std::vector<double> freq, int c,bool instantPlot)
{
	JSL::gnuplot gp;
	gp.WindowSize(1000,800);
	basicPlot(gp,d,c);
	
	TransitionPlot(gp,idx,freq,d.Chromosomes[c].maxIdx);

	gp.SetTitle(d.Chromosomes[c].Name);
	std::string outname = settings.OutputName + "_" + d.Chromosomes[c].Name + ".png";
	gp.SetXLabel("Chromosome Index (bp)");
	gp.SetYLabel("Coverage");
	gp.SetLegend(true);
	int maxVal = std::min((int)(d.Mean + 3*d.Deviation),d.Chromosomes[c].maxK);
	auto it = max_element(std::begin(freq), std::end(freq)); 
	maxVal = std::max(maxVal,(int)(*it*1.2));
	gp.SetYRange(0,maxVal);
	if (!instantPlot)
	{
		gp.SetOutput(outname);
		gp.SetTerminal("png");
	}
	gp.Show();
}

inline void OutputPlot(Data & d, Settings & settings)
{
	//write it to run from .tree file so can be run independently of the code which generates Transition objects -- makes it a bit roundabout, but will be more useful in the long run!
	Log("\tBeginning Plot Routine\n");
	std::string reconstructFile =settings.OutputName + ".tree";
	double nu;
	double sigma;	
	double gamma;
	double L;
	bool inMainBody = false;
	std::vector<std::vector<chr_int>> idx;
	std::vector<std::vector<double>> freq;
	std::vector<int> trueChromIdx;
	std::string currentChrom = "";
	int chromID=-1;
	bool allowed = false;
	forLineVectorIn(reconstructFile,' ',
		if (inMainBody)
		{
			auto chrom = FILE_LINE_VECTOR[0];
			if (chrom != currentChrom)
			{
				currentChrom = chrom;
				
				if (currentChrom == settings.TargetChromosome || settings.TargetChromosome == "all")
				{
					++chromID;
					idx.resize(chromID+1);
					freq.resize(chromID+1);

					int trueChrom = -1;
					for (int c = 0; c < d.Chromosomes.size(); ++c)
					{
						if (d.Chromosomes[c].Name == currentChrom)
						{
							trueChrom = c;
							break;
						}
					}
					if (trueChrom == -1)
					{
						std::cout << "I could not locate " << currentChrom << " in my raw datafile -- this indicates that the processed data and the raw data are from different samples!" << std::endl;
						exit(2);
					}
					trueChromIdx.push_back(trueChrom);
					allowed = true;
				}
				else
				{
					allowed = false;
				}
			}

			if (allowed)
			{
				idx[chromID].push_back(std::stoi(FILE_LINE_VECTOR[1]));
				freq[chromID].push_back(nu*std::stoi(FILE_LINE_VECTOR[2]));
			}
			
		}
		else
		{
			if (FILE_LINE_VECTOR.size() < 2)
			{
				inMainBody = true;
			}
			else
			{
				if (FILE_LINE_VECTOR[0] == "Nu")
				{
					nu = std::stod(FILE_LINE_VECTOR[1]);
				}
				if (FILE_LINE_VECTOR[0] == "Sigma")
				{
					sigma = std::stod(FILE_LINE_VECTOR[1]);
				}
				if (FILE_LINE_VECTOR[0] == "Gamma")
				{
					gamma = std::stod(FILE_LINE_VECTOR[1]);
				}
				if (FILE_LINE_VECTOR[0] == "L")
				{
					L = std::stod(FILE_LINE_VECTOR[1]);
				}
			}
		}
	);

	for (int i = 0; i < idx.size(); ++i)
	{
		Log("\t\tPlotting " << d.Chromosomes[trueChromIdx[i]].Name << "\n");
		TransitionFrame(d,settings,idx[i],freq[i],trueChromIdx[i],settings.PlotOnly);
	}
	Log("\tPlotting Completed\n");
}