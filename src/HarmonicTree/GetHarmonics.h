#pragma once
#include "../DataFrame/data.h"
#include "../settings.h"
#include <vector>
#include <chrono>
#include "HarmonicNetwork.h"
#include "../Utility/Distributor.h"

Path GetHarmonics(Data & d, Settings & s,JSL::gnuplot & gp);