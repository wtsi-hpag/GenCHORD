#include "plotting.h"

 void basicPlot(JSL::gnuplot & gp,Data & d,int chrom)
{
	gp.Plot(d.Chromosomes[chrom].Idx,d.Chromosomes[chrom].Counts,JSL::LineProperties::Legend("Data"),JSL::LineProperties::Colour("black"));
	
}

 void TransitionPlot(JSL::gnuplot & gp,std::vector<chr_int> idx, std::vector<double> freq,double nuCorrect,std::string leg)
{
	std::vector<chr_int> remap_idx;
	std::vector<double> remap_freq;
	double prev = freq[0]/nuCorrect;
	for (int i = 0; i < idx.size(); ++i)
	{
		remap_idx.push_back(idx[i]);
		remap_idx.push_back(idx[i]);
		remap_freq.push_back(prev);
		prev = freq[i]/nuCorrect;
		remap_freq.push_back(prev);
	}

	gp.Plot(remap_idx,remap_freq,JSL::LineProperties::Legend(leg),JSL::LineProperties::PenSize(2));
}

TreeMeta::TreeMeta(Settings & settings)
{
	auto q = JSL::split(settings.DataFile,'.');
	std::string reconstructFile = q[0] + ".tree";
	settings.OutputName = q[0];


	bool inMainBody = false;
	std::vector<std::vector<chr_int>> idx;
	std::vector<std::vector<double>> freq;
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
					Idx.resize(chromID+1);
					Frequency.resize(chromID+1);
					Harmonic.resize(chromID+1);
					Name.push_back(currentChrom);
					allowed = true;
				}
				else
				{
					allowed = false;
				}
			}

			if (allowed)
			{
				Idx[chromID].push_back(std::stoi(FILE_LINE_VECTOR[1]));
				Frequency[chromID].push_back(Nu*std::stoi(FILE_LINE_VECTOR[2]));
				Harmonic[chromID].push_back(std::stoi(FILE_LINE_VECTOR[2]));
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
				if (FILE_LINE_VECTOR[0] == "Data")
				{
					DataFile = FILE_LINE_VECTOR[1];
				}
				if (FILE_LINE_VECTOR[0] == "Nu")
				{
					Nu = std::stod(FILE_LINE_VECTOR[1]);
				}
				if (FILE_LINE_VECTOR[0] == "Sigma")
				{
					Sigma = std::stod(FILE_LINE_VECTOR[1]);
				}
				if (FILE_LINE_VECTOR[0] == "Gamma")
				{
					Gamma = std::stod(FILE_LINE_VECTOR[1]);
				}
				if (FILE_LINE_VECTOR[0] == "L")
				{
					L = std::stoi(FILE_LINE_VECTOR[1]);
				}
			}
		}
	);
}

 void TransitionFrame(Data & d, Settings & settings, std::vector<chr_int> idx, std::vector<double> freq, int c,bool instantPlot)
{
	JSL::gnuplot gp;
	gp.WindowSize(1000,800);
	basicPlot(gp,d,c);
	
	TransitionPlot(gp,idx,freq,1,"Treefit");

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

 void OutputPlot(Settings & settings)
{
	//write it to run from .tree file so can be run independently of the code which generates Transition objects -- makes it a bit roundabout, but will be more useful in the long run!
	Log("\tBeginning Plot Routine\n");
	std::string reconstructFile =settings.OutputName + ".tree";
	TreeMeta T(settings);
	Data d(T.DataFile,settings.DataThinning,settings.TargetChromosome,settings.MemorySmoothing);	
	std::vector<int> nameRemapper;
	for (int i = 0; i < T.Name.size(); ++i)
	{
		int remap = -1;
		for (int j = 0; j < d.Chromosomes.size(); ++j)
		{
			if (d.Chromosomes[j].Name == T.Name[i])
			{
				remap = j;
			}
		}
		if (remap == -1)
		{
			std::cout << " I seem to have encountered an error mapping Tree files to the data" << std::endl;
			exit(5);
		}
		else
		{
			nameRemapper.push_back(remap);
		}
	}

	for (int i = 0; i < T.Idx.size(); ++i)
	{
		Log("\t\tPlotting " << T.Name[i] << "  " << nameRemapper[i] << "\n");
		TransitionFrame(d,settings,T.Idx[i],T.Frequency[i],nameRemapper[i],settings.PlotOnly);
	}
	Log("\tPlotting Completed\n");
}

void TreePlot(JSL::gnuplot & gp,TreeMeta T,int j,std::vector<int> & maxQ)
{

	// for (int T.Idx)

	auto d = std::max_element(T.Frequency[j].begin(),T.Frequency[j].end());
	int mQ = round(*d/T.Nu);

	if (j >= maxQ.size())
	{
		maxQ.resize(j+1,0);
	}
	if (mQ > maxQ[j])
	{
		maxQ[j] = mQ;

	}
	TransitionPlot(gp,T.Idx[j],T.Frequency[j],T.Nu,T.DataFile);


}

void ComparisonPlots(Settings & settings)
{
	auto dirSplit = JSL::split(settings.OutputName,'/');
	std::string frontLoad = "";
	if (dirSplit.size() > 1)
	{
		for (int i =0; i < dirSplit.size(); ++i)
		{
			frontLoad += dirSplit[i] + "/";
		}
		JSL::mkdir(frontLoad);
	}


	forLineVectorIn(settings.ComparePlot,' ',
		std::string title = "";
		auto first = FILE_LINE_VECTOR[0];
		int nid = 0;
		if (first.substr(0,1)=="#")
		{
			nid = 1;
			title = first.substr(1,std::string::npos) + "-";
		}
		settings.DataFile = FILE_LINE_VECTOR[nid];
		std::vector<JSL::gnuplot> gps;
		std::vector<std::string> names;
		std::vector<int> maxQ(FILE_LINE_VECTOR.size()+2,0);
		for (int i = nid; i < FILE_LINE_VECTOR.size(); ++i)
		{
			auto file = FILE_LINE_VECTOR[i];
			settings.DataFile = file;
			TreeMeta T(settings);
			if (T.Name.size() > gps.size())
			{
				gps.resize(T.Name.size());
			}
			
			names = T.Name;
			for (int j = 0; j < T.Name.size(); ++j)
			{
				// gps[j].Plot(T,j);
				// TreePlot(gps[j],T,j,maxQ);
				gps[j].SetTitle(title + T.Name[j]);
				
			}
		}


		for (int i = 0; i < gps.size(); ++i)
		{
			gps[i].SetXLabel("Chromosome Index (bp)");
			gps[i].SetYLabel("Harmonic");
			std::string n = frontLoad + title+names[i] +".png";
			gps[i].SetOutput(n);
			gps[i].SetTerminal("png");
			gps[i].SetLegend(true);
			gps[i].SetYRange(0,maxQ[i]+0.5);
			// std::cout << "Saving to " << n << "  " << maxQ[i] << "  " << maxQ.size() << "  " << i << std::endl;
			
			gps[i].Show();
		}

	);

}