#pragma once
#include "settings.h"
#include "Transitions.h"
#include "data.h"
#include "NegativeBinomial.h"

class ErroredBinomial;

extern int Ploidy;
extern double logPloidyPrior;
struct qReturn
{
	chr_int idx;
	int Q;
	double Score;
	qReturn(int q, double s): Q(q), Score(s){};
};

qReturn bestQ(int chrom,int iStart,long long int maxIdx, double currentQ, double alpha, ErroredBinomial * EB, Data & d, int accelerator);

qReturn EdgeFinder(int chrom,int iStart, int L, int originalQ, Data & d, ErroredBinomial * EB, double alpha, int accelerator);

void ChromosomeAssign(Transitions * output, int chrom, Data & d,double alpha, ErroredBinomial * EB, int L, int accelerator,double nu);