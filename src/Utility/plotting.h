#pragma once
#include "../../libs/JSL/JSL.h"
#include "../data.h"
#include "../Transitions.h"
#include "../settings.h"
#include "../HarmonicTree/Path.h"
extern int nPlot;
void basicPlot(JSL::gnuplot & gp,Data & d,int chrom);

void TransitionPlot(JSL::gnuplot & gp,const Data & d, Path path, std::string legend);
