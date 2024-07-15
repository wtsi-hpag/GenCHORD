#include "plotting.h"

 void basicPlot(JSL::gnuplot & gp,Data & d,int chrom)
{
	gp.Plot(d.Chromosomes[chrom].Idx,d.Chromosomes[chrom].SmoothCounts,JSL::LineProperties::Legend("Data"),JSL::LineProperties::Colour("black"));
	
}

void TransitionPlot(JSL::gnuplot & gp,const ChromosomeCoverage & chrom, Path & path, std::string legend)
{
	std::vector<int> plotIdx = {path.Route[0].Index};
	int q = path.Route[0].Value;
	double qPrime = (1.0 - path.Eta[0])*q + path.Eta[0]*path.Dc;
	double prev = q*path.Nu;
	double contprev = qPrime*path.Nu;
	std::vector<double> plotVals = {prev};
	std::vector<double> contVals = {contprev};
	for (int i =1; i < path.Route.size(); ++i)
	{

		plotIdx.push_back(path.Route[i].Index);
		plotVals.push_back(prev);
		contVals.push_back(contprev);
		int q = path.Route[i].Value;
		double qPrime = (1.0 - path.Eta[q])*q + path.Eta[q]*path.Dc;
		contprev = qPrime*path.Nu;
		prev = q * path.Nu;
		plotIdx.push_back(path.Route[i].Index);
		plotVals.push_back(prev);
		contVals.push_back(contprev);
		// plotIdx.push_back(path.Route[i-1].Index);
		// plotIdx.push_back(path.Route[i].Index);
		//  prev = path.Route[i].Value*path.Nu;
		//  plotVals.push_back(prev);
		//  plotVals.push_back(prev);
	}
	plotIdx.push_back(chrom.maxIdx);
	plotVals.push_back(prev);
	contVals.push_back(contprev);
	gp.Plot(plotIdx,plotVals,JSL::LineProperties::Legend(legend),JSL::LineProperties::PenSize(6));
	gp.Plot(plotIdx,contVals,JSL::LineProperties::Legend("Including Contamination"),JSL::LineProperties::PenSize(3));
	gp.SetXLabel("Chromosome Index");
	gp.SetYLabel("Coverage");
}