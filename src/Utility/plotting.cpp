#include "plotting.h"

 void basicPlot(JSL::gnuplot & gp,Data & d,int chrom)
{
	gp.Plot(d.Chromosomes[chrom].Idx,d.Chromosomes[chrom].Counts,JSL::LineProperties::Legend("Data"),JSL::LineProperties::Colour("black"));
	
}

void TransitionPlot(JSL::gnuplot & gp,const ChromosomeCoverage & chrom, Path & path, std::string legend)
{
	std::vector<int> plotIdx = {path.Route[0].Index};
	double prev = path.Route[0].Value * path.Nu;
	std::vector<double> plotVals = {prev};
	for (int i =1; i < path.Route.size(); ++i)
	{
		plotIdx.push_back(path.Route[i].Index);
		plotVals.push_back(prev);
		prev = path.Route[i].Value*path.Nu;
		plotIdx.push_back(path.Route[i].Index);
		plotVals.push_back(prev);
	}
	plotIdx.push_back(chrom.maxIdx);
	plotVals.push_back(prev);

	gp.Plot(plotIdx,plotVals,JSL::LineProperties::Legend(legend),JSL::LineProperties::PenSize(3));
}