#pragma once

#include "mc.h"

void gibbsVolumeMove(MonteCarlo& systemA, MonteCarlo& systemB);
void gibbsSwapMove(MonteCarlo& systemA, MonteCarlo& systemB);
void runGibbsEnsemble(MonteCarlo& systemA, MonteCarlo& systemB);