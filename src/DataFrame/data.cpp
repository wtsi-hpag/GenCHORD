#include "data.h"


Data::Data()
{
	Loaded = false;
}

Data::Data(Settings & settings)
{
	PrepareForLoading();
	if (JSL::PipedInputFound())
	{
		Log("Loading Data in from pipe\n");
		forLineVectorInPipedInput(settings.DataFileDelimiter,
			ParseLine(PIPE_LINE_VECTOR,settings);
		)

		if (File.ChromActive)
		{
			Chromosomes.push_back(File.Chromosome);
		}
	}
	else
	{
		if (settings.DataFile != "_no_file_provided_")
		{
			Log("Loading Data in from " << settings.DataFile << "\n");
			forLineVectorIn(settings.DataFile,settings.DataFileDelimiter,
				ParseLine(FILE_LINE_VECTOR,settings);
			);
			if (File.ChromActive)
			{
				Chromosomes.push_back(File.Chromosome);
			}
		}
		else
		{
			std::cout << "\nERROR No Data detected. Please provide a valid data input\n" << std::endl;
			exit(1);
		}
	}

	Mean = File.CoverageSum/File.LoadedData;
	Deviation = sqrt(File.CoverageSquareSum/ File.LoadedData - Mean * Mean);

	//do some corrections
	int LudicrousValue = Mean + 10*Deviation;
	int truncated = 0;
	for (int i = 0; i < Chromosomes.size(); ++i)
	{
		truncated += Chromosomes[i].RemoveSpuriousData(LudicrousValue);
		File.MaxCoverage = std::max(File.MaxCoverage,Chromosomes[i].maxK);
	}
	maxK = std::min(LudicrousValue,File.MaxCoverage);
	Log("\tAsserting that coverage above " << LudicrousValue << " is spurious. " << truncated << " Datapoints were affected\n")
	Log("\tSpoofed " << File.SpoofedData << " missing entries as gaps\n\tA total of " << File.LoadedData << " datapoints loaded\n\tData has global coverage  "<< Mean << "±" << Deviation  << "\n\tMaximum coverage value is " << File.MaxCoverage<< std::endl;);


	Loaded = true;
	std::cout << "\t" << File.LoadedData << " points loaded in total" << std::endl;
}

void Data::PrepareForLoading()
{
	File.CurrentLine = 0;
	File.CurrentChromosome = 0;
	File.ChromActive = false;
	File.LoadedData = 0;
	File.SpoofedData = 0;
	File.MaxCoverage = 0;
	File.CoverageSquareSum =0;
	File.CoverageSum = 0;
}	


void Data::AddData(int index, int coverage, const Settings & settings,bool spoofed)
{
	File.SmoothingAverage = settings.MemorySmoothing * File.SmoothingAverage + (1.0 - settings.MemorySmoothing) * coverage;
	
	if (index % settings.DataThinning == 0)
	{
		int k = round(File.SmoothingAverage);
		
		File.Chromosome.Add(index,k);
		File.CoverageSum += File.SmoothingAverage;
		File.CoverageSquareSum += pow(File.SmoothingAverage,2);
		++File.LoadedData;
		if (spoofed)
		{
			++File.SpoofedData;
		}
	}
	
}

void Data::ParseLine(const std::vector<std::string> & line, const Settings & settings)
{
	//check which chromosome this line is associated with
	if (File.CurrentLine == 0)
	{
		File.ChromSignifier = line[0];
		File.PreviousIndex = std::stoi(line[1]);

		if (settings.AllChromosomes || File.ChromSignifier == settings.TargetChromosome)
		{
			// Log("\tFound Chromosome " << File.CurrentChromosome+1 << " on line " << File.CurrentLine << "\n");
			File.ChromActive = true;
			File.Chromosome = ChromosomeCoverage(File.ChromSignifier);
			File.SmoothingAverage = std::stoi(line[2]);
		}
	}
	else
	{
		if (line[0] != File.ChromSignifier)
		{
			File.ChromSignifier = line[0];
			++File.CurrentChromosome;
			File.PreviousIndex = std::stoi(line[1]) - File.GapSize; //spoof in the jumps in index that occur when moving to a new chromosome!
			File.SmoothingAverage = std::stoi(line[2]);
			if (File.ChromActive)
			{
				Chromosomes.push_back(File.Chromosome);
			}

			if (settings.AllChromosomes || File.ChromSignifier == settings.TargetChromosome)
			{
				// Log("\tFound Chromosome " << File.CurrentChromosome+1 << " on line " << File.CurrentLine<<"\n");
				File.ChromActive = true;
				File.Chromosome = ChromosomeCoverage(File.ChromSignifier);
			}
			else
			{
				File.ChromActive= false;
			}
		}
	}

	//detect gap size -- assumes that line 2 and line 1 have no omitted datapoints
	if (File.CurrentLine == 1)
	{
		File.GapSize = std::stoi(line[1]) - File.PreviousIndex;
	}

	//start loading up the data
	if (File.ChromActive)
	{
		chr_int index = std::stoi(line[1]);
		int coverage = std::stoi(line[2]);
		if (coverage > 2000)
		{
			coverage = File.SmoothingAverage;
		}

		
		//detect gaps in the file, which occur when many consectutive zeros
		if (File.CurrentLine > 1)
		{
			if (index - File.PreviousIndex < File.GapSize)
			{
				File.GapSize = index - File.PreviousIndex;
			}
			while (File.GapSize + File.PreviousIndex < index)
			{
				File.PreviousIndex += File.GapSize;

				AddData(File.PreviousIndex,0,settings,true);
			
			}
		}
		AddData(index,coverage,settings,false);

		File.PreviousIndex = index;

	}

	++File.CurrentLine;
}
