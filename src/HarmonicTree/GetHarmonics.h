#pragma once
#include "../DataFrame/data.h"
#include "../settings.h"
#include <vector>
#include <chrono>
#include "HarmonicNetwork.h"
#include "../Utility/Distributor.h"
#include "../Probability/ProbabilityModel.h"

std::vector<Path> GetHarmonics(Data & d, Settings & s,std::vector<ProbabilityModel> & prob);