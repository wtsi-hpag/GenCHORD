#pragma once
#include "../data.h"
#include "../../libs/JSL/JSL.h"
#include "../Probability/ErroredBinomial.h"
#include "../Utility/plotting.h"
#include "../settings.h"
#include "../Transitions.h"
#include "HarmonicFunctions.h"
void HarmonicFit(Data & d, const Settings & settings);