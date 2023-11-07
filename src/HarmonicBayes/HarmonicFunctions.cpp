#include "HarmonicFunctions.h"

int Ploidy = 2;
double logPloidyPrior = log(0.5);
qReturn bestQ(int chrom,int iStart,long long int maxIdx, double currentQ, double alpha, ErroredBinomial * EB, Data & d, int accelerator)
{
	double bestQ = -1;
	double bestVal = -99999999999;
	chr_int maxI = d.Chromosomes[chrom].Idx.size();
	for (int q = 0; q < EB->qMax; ++q)
	{
		double logp = 0;
		if (q != currentQ)
		{
			logp += log(alpha);
		}
		
		int i = iStart;
		chr_int gidx = d.Chromosomes[chrom].Idx[i];
		while (gidx < maxIdx)
		{
			if (q != Ploidy)
			{
				logp += logPloidyPrior;
			}
			logp += EB->GetProb(q,d.Chromosomes[chrom].Counts[i]);
			i+=accelerator;
			if ( i < maxI)
			{
				gidx = d.Chromosomes[chrom].Idx[i];
			}
			else
			{
				gidx =maxIdx + 10;
			}
		}

		if (q == 0 || logp > bestVal)
		{
			bestQ = q;
			bestVal = logp;
		}
	}
	if (bestQ == currentQ)
	{
		bestVal = EB->GetProb(currentQ,d.Chromosomes[chrom].Counts[iStart]); 
		if (currentQ != Ploidy)
		{
			bestVal += logPloidyPrior;
		}
		// hack such that if no transition, returns only the score of the ith base, rather than the [i,i+L] domain -- allows accumulation
	}
	return qReturn(bestQ,bestVal);
}


qReturn EdgeFinder(int chrom,int iStart, int L, int originalQ, Data & d, ErroredBinomial * EB, double alpha, int accelerator)
{
	// std::cout << "Edge finding" << "  " << iStart << "  " << d.Chromosomes[chrom].Idx.size() << std::endl;
	int i = iStart;
	double tail = 0;//EB.GetProb(originalQ,d.counts[i]);
	int offset =0;
	chr_int idx = d.Chromosomes[chrom].Idx[i+offset];
	chr_int stop = idx + L; //search up to half the length of the gap -- but any further and the gap is just being mollycoddled and your prior should have caught it
	chr_int search_Stop = idx + L/2;
	qReturn out = bestQ(chrom,i+offset,stop,originalQ,alpha,EB,d,accelerator);
	out.idx = idx;
	int maxI = d.Chromosomes[chrom].Idx.size();
	double ploidyError = 0;
	if (originalQ != Ploidy)
	{
		ploidyError = logPloidyPrior;
	}
	// std::cout << "\t" << originalQ << "-> " << out.Q << out.Score << std::endl;
	while(idx < search_Stop && idx < d.Chromosomes[chrom].maxIdx)
	{
		auto refinedInfer = bestQ(chrom,i+offset,stop,originalQ,alpha,EB,d,accelerator);
		refinedInfer.Score += tail;

		double t =  EB->GetProb(originalQ,d.Chromosomes[chrom].Counts[i+offset]) + ploidyError;
		tail += t;

		if (refinedInfer.Score > out.Score && refinedInfer.Q != originalQ) 
		{
			out = refinedInfer;
			out.idx = idx;
		}
		offset+=accelerator;
		if (i+offset < maxI)
		{
			idx = d.Chromosomes[chrom].Idx[i+offset];
		}
		else
		{
			idx = search_Stop +4;
			// std::cout << "oh no" << std::endl;
		}
	}
	// exit(5);
	// std::cout << "Found edge" << std::endl;
	return out;
}

std::mutex qmtx;
void ChromosomeAssign(Transitions * output, int chrom, Data & d,double alpha, ErroredBinomial * EB, int L, int accelerator,double nu)
{
	chr_int maxIndex = d.Chromosomes[chrom].maxIdx;
		
	int q = -1;
	chr_int index = 0;
	int i = 0;
	double s = 0;
	while (i < d.Chromosomes[chrom].Idx.size())
	{
		index = d.Chromosomes[chrom].Idx[i];
		auto infer = bestQ(chrom,i,index+L,q,alpha,EB,d,accelerator);
		s += infer.Score;
		
		if (infer.Q != q)
		{
			if (i > 0)
			{
				// std::cout << q << "-> " << infer.Q << infer.Score << std::endl;
				auto refinedInfer = EdgeFinder(chrom,i,L,q,d,EB,alpha,accelerator);

				qmtx.lock();
				output->Add(chrom,refinedInfer.idx,nu,EB->Sigma,refinedInfer.Q);
				qmtx.unlock();

				q = refinedInfer.Q;
			}
			else
			{
				qmtx.lock();
				output->Add(chrom,index,nu,EB->Sigma,infer.Q);
				qmtx.unlock();
				q = infer.Q;
			}
			chr_int stop = d.Chromosomes[chrom].Idx[i] + L;
			while (index < stop)
			{
				i+=accelerator;
				if (i < d.Chromosomes[chrom].Idx.size())
				{
					index = d.Chromosomes[chrom].Idx[i];
				}
				else
				{
					index = stop + 10;
					i = d.Chromosomes[chrom].Idx.size()+5;
				}
			}
		}
		else
		{
			i+=accelerator;
		}

	}
	qmtx.lock();
	output->Score += s;
	qmtx.unlock();
}