#pragma once

#include "../DataFrame/data.h"
#include "../Probability/ProbabilityModel.h"




std::vector<ProbabilityModel> NormaliseModel(ProbabilityModel & p, const Data & data, const Settings & settings, int kMax, int Qmax);
