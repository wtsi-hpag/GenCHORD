#include "Distributor.h"

std::mutex mtx;

struct ShuffReturn
{
	std::vector<std::vector<int>> Body;
	std::vector<int> Header;

	double variance(Data & d)
	{
		double sqSum = 0;
		double sum = 0;
		for (int i = 0; i < Body.size(); ++i)
		{
			double L = 0;
			for (int j = 0; j < Body[i].size(); ++j)
			{
				L+= d.Chromosomes[Body[i][j]].maxIdx;
			}
			sum += L;
			sqSum += L*L; 
		// std::cout << "\tWorker " << i << " has vector " << JSL::Vector(Body[i]) << " L: " << L << std::endl;
		}
		double L = 0;
		for (int j = 0; j < Header.size(); ++j)
		{
			L += d.Chromosomes[Header[j]].maxIdx;
		}
		// std::cout << "\tMain thread has vector " << JSL::Vector(Header) << " L: " << L << std::endl;
		sum += L;
		sqSum += L*L; 
		int N = (Body.size()+1);
		double mean = (double)sum/N;
		double sqMean = (double)sqSum/N;

		double std = sqrt(sqMean - mean*mean)/mean;
		return std;
	}
};

ShuffReturn OptimizeDistrubution(Data & d, int nWork)
{
	std::vector<int> cs(d.Chromosomes.size());
	std::iota(cs.begin(),cs.end(),0);
	 std::random_device rd;
    std::mt19937 g(rd());
	srand(time(NULL));
	ShuffReturn out;
	double var = 9999999;
	int attempts = 50;
	for(int l = 0; l < attempts; ++l)
	{
		// std::shuffle(cs.begin(),cs.end(),g);
		ShuffReturn attempt;
		attempt.Body.resize(nWork);

		int insert = 0;
		for (int i = 0; i < cs.size(); ++i)
		{
			if (l > attempts/2*10)
			{
				insert = rand() % (nWork + 1); //really go hog wild
			}		
			if (insert < nWork)
			{
				attempt.Body[insert].push_back(cs[i]);
				++insert;
			}
			else
			{
				attempt.Header.push_back(cs[i]);
				insert = 0;
			}
		}
		double attemptVar = attempt.variance(d);
		if (attemptVar < var)
		{
			out = attempt;
			var = attemptVar;
		}
	}
	if (nWork > 0)
		Log("\t\tDistributed a core-deviance of " << var << std::endl;)
	return out;
}


Distributor::Distributor(int nWorkers, Data & d) : DataCopy(d)
{
	if (nWorkers > 0)
		Log("\tSetting up asynchronous pool"<<std::endl;);
	WorkerCount = nWorkers;
	
	ChromAssignments.resize(nWorkers);
	
	auto S = OptimizeDistrubution(d,nWorkers);
	ChromAssignments = S.Body;
	MainChromAssigment = S.Header;

	


	Register = std::vector<int>(nWorkers,-2);
	ebStart.resize(WorkerCount);
	ebEnd.resize(WorkerCount);
	for (int i = 0; i < WorkerCount; ++i)
	{
		Threads.push_back(std::thread(&Distributor::Worker,this,i));
		// Threads[i].detach();
	}
	Gather();
	Log("\t\tAll workers initialised & ready to work\n");
}

void Distributor::UpdateEB(ErroredBinomial * eb)
{
	int delta = eb->Resolution/(WorkerCount+1);
	EB = eb;
	for (int i = 0; i < WorkerCount; ++i)
	{
		ebStart[i] = i*delta;
		ebEnd[i] = (i+1)*delta;
	}
}
void Distributor::SpoolDown(int id)
{
	const std::lock_guard<std::mutex> lock(mtx);
	Register[id] = 0;

}

int Distributor::CheckRegister(int id)
{
	const std::lock_guard<std::mutex> lock(mtx);
	int r = Register[id];
	return r;
}
void Distributor::SetRegister(int value, int id)
{
	const std::lock_guard<std::mutex> lock(mtx);
	Register[id] = value;
}
void Distributor::Worker(int id)
{
	mtx.lock();
	// Log("\t\tWorker " << id << " initialised " << std::endl;);
	Register[id] = 0;
	mtx.unlock();
	bool WorkerActive = true;
	while (WorkerActive)
	{
		
		int q = CheckRegister(id);

		switch (q)
		{
			case -1:
			{
				WorkerActive = false;
				SpoolDown(id);
				// return;
				break;
			}
			case 1:
			{
				EB_Task(id);
				break;
			}
			case 2:
			{
				Assignment_Task(id);
				break;
			}
			default:
			{
				break;
			}
		}

	}

}

void Distributor::EB_Task(int id)
{
	EB->PopulateChunk(ebStart[id],ebEnd[id]);
	SetRegister(0,id);
}

void Distributor::Assignment_Task(int id)
{
	for (int chrom : ChromAssignments[id])
	{
		ChromosomeAssign(Output,chrom,DataCopy,alpha,EB,L,accelerator,nu);
	}
	SetRegister(0,id);
}
void Distributor::UpdateAssigner(Transitions * output, double a, int l, int acc, double freq)
{
	Output = output;
	alpha = a;
	L = l;
	accelerator = acc;
	nu = freq;
}
void Distributor::Signal(int value)
{
	for (int i = 0; i < Register.size(); ++i)
	{
		int q = CheckRegister(i);
		if (q == 0)
		{
			SetRegister(value,i);
		}
		else
		{
			std::cout << "Error! Tried to dispatch worker " << i << " on task " << value << ", but previous task (" << q << ") not completed" << std::endl;
			exit(5);
		}
	}
}
void Distributor::Gather()
{
	int completed = 0;
	while (completed != WorkerCount)
	{
		completed = 0;
	
		for (int j = 0; j < WorkerCount; ++j)
		{
			int q = CheckRegister(j);
			if (q == 0)
			{
				++completed;
			}
		}
	}
}
Distributor::~Distributor()
{
	
	Signal(-1);
	Gather();
	for (int i = 0; i < WorkerCount; ++i)
	{
		if (Threads[i].joinable())
		{
			Threads[i].join();
		}
	}
	Threads.resize(0);
	Log("Async Pool terminated\n");
}
