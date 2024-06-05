#pragma once

#include "../DataFrame/data.h"
#include "../Probability/ProbabilityModel.h"




void NormaliseModel(ProbabilityModel & p, const Data & data, const Settings & settings, int kMax, int Qmax);
