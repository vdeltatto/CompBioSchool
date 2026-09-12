
#include <cmath>
#include <random>
#include <vector>

#include "mc.h"

bool parallelTemperingSwap(MonteCarlo& a, MonteCarlo& b)
{
  double acc = std::exp(-(b.beta - a.beta) * (b.runningEnergies.energy - a.runningEnergies.energy));

  if (uniform() < acc)
  {
    std::swap(a.temperature, b.temperature);
    std::swap(a.beta, b.beta);
    return true;
  }
  return false;
}