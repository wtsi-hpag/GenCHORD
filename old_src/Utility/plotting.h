#pragma once
#include "../../libs/JSL/JSL.h"
#include "../DataFrame/data.h"
#include "../Transitions.h"
#include "../settings.h"
#include "../HarmonicTree/Path.h"
extern int nPlot;
void basicPlot(JSL::gnuplot & gp,Data & d,int chrom);

void TransitionPlot(JSL::gnuplot & gp,const ChromosomeCoverage & chrom, Path & path, std::string legend);
