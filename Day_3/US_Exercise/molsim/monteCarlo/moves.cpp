#include <cmath>
#include <vector>

#include "mc.h"

void MonteCarlo::updateMoveStatistics()
{
  translationAcceptance = translationAttempted ? translationAccepted / translationAttempted : 0.0;
  volumeAcceptance = volumeAttempted ? volumeAccepted / volumeAttempted : 0.0;
  insertionAcceptance = insertionAttempted ? insertionAccepted / insertionAttempted : 0.0;
  deletionAcceptance = deletionAttempted ? deletionAccepted / deletionAttempted : 0.0;
}

void MonteCarlo::translationMove()
{
  translationAttempted++;
  // Generate random displacement
  int particleIdx = static_cast<int>(uniform() * numberOfParticles);
  double3 displacement = maxDisplacement * (double3(uniform(), uniform(), uniform()) - 0.5);
  double3 trialPosition = positions[particleIdx] + displacement;

  // Calculate old and new energy and virial
  EnergyVirial oldEnergyVirial = particleEnergyVirial(positions, positions[particleIdx], particleIdx, boxSize, cutOff);
  EnergyVirial newEnergyVirial = particleEnergyVirial(positions, trialPosition, particleIdx, boxSize, cutOff);

  // Accept or reject the move based on Metropolis criterion
  if (uniform() < std::exp(-beta * (newEnergyVirial.energy - oldEnergyVirial.energy)))
  {
    translationAccepted++;
    positions[particleIdx] = trialPosition;
    runningEnergyVirial += (newEnergyVirial - oldEnergyVirial);
  }
}

void MonteCarlo::optimizeMaxDisplacement()
{
  // Adjust max displacement to maintain optimal acceptance ratio
  if (translationAttempted > 100)
  {
    double acceptance = translationAccepted / translationAttempted;
    double scaling = std::clamp(2.0 * acceptance, 0.5, 1.5);
    maxDisplacement = std::clamp(scaling * maxDisplacement, 0.0001, 0.49 * boxSize);
    translationAccepted = 0;
    translationAttempted = 0;
  }
}

void MonteCarlo::volumeMove()
{
  volumeAttempted++;

  double volumeChange = (uniform() - 0.5) * maxVolumeChange;

  double newVolume = volume + volumeChange;
  if (newVolume < 0.0)
  {
    return;
  }

  double newBoxSize = std::cbrt(newVolume);
  double scale = newBoxSize / boxSize;

  std::vector<double3> trialPositions(positions);
  for (int i = 0; i < numberOfParticles; i++)
  {
    trialPositions[i] *= scale;
  }
  EnergyVirial newEnergyVirial = systemEnergyVirial(trialPositions, newBoxSize, cutOff, cutOffPrefactor);

  // start refactor
  if (uniform() < 0.0)
  // end refactor
  {
    volumeAccepted++;
    positions = trialPositions;
    boxSize = newBoxSize;
    volume = boxSize * boxSize * boxSize;
    density = numberOfParticles / volume;
    runningEnergyVirial = newEnergyVirial;
  }
}

void MonteCarlo::optimizeVolumeChange()
{
  if (volumeAttempted > 100)
  {
    double acceptance = volumeAccepted / volumeAttempted;
    double scaling = std::clamp(2.0 * acceptance, 0.5, 1.5);
    maxVolumeChange = std::clamp(scaling * maxVolumeChange, 0.001 * volume, 0.5 * volume);
    volumeAccepted = 0;
    volumeAttempted = 0;
  }
}

void MonteCarlo::swapMove()
{
  if (uniform() < 0.5)
  {
    insertionAttempted++;
    // trial insertion of new particle
    double3 newParticle = boxSize * double3(uniform(), uniform(), uniform());
    EnergyVirial diffEnergyVirial = particleEnergyVirial(positions, newParticle, numberOfParticles, boxSize, cutOff);

    // add increase tail energy (N+1)**2 - N**2
    diffEnergyVirial += (cutOffEnergyVirial(numberOfParticles + 1, boxSize, cutOffPrefactor) -
                         cutOffEnergyVirial(numberOfParticles, boxSize, cutOffPrefactor));

    double Pacc = (pressure * beta * volume / static_cast<double>(numberOfParticles + 1.0)) *
                  std::exp(-beta * diffEnergyVirial.energy);

    if (uniform() < Pacc)
    {
      // accept
      insertionAccepted++;
      runningEnergyVirial += diffEnergyVirial;
      numberOfParticles++;
      density = numberOfParticles / volume;
      positions.push_back(newParticle);
    }
  }
  else
  {
    if (numberOfParticles == 0)
    {
      // can not delete if there are no particles
      return;
    }

    deletionAttempted++;
    // trial deletion
    int selectedParticle = static_cast<int>(uniform() * numberOfParticles);
    EnergyVirial diffEnergyVirial =
        particleEnergyVirial(positions, positions[selectedParticle], selectedParticle, boxSize, cutOff);
    diffEnergyVirial *= -1.0;

    // removed tail energy
    diffEnergyVirial += (cutOffEnergyVirial(numberOfParticles - 1, boxSize, cutOffPrefactor) -
                         cutOffEnergyVirial(numberOfParticles, boxSize, cutOffPrefactor));

    double Pacc = (numberOfParticles / (pressure * beta * volume)) * std::exp(-beta * diffEnergyVirial.energy);

    if (uniform() < Pacc)
    {
      deletionAccepted++;
      runningEnergyVirial += diffEnergyVirial;
      numberOfParticles--;
      density = numberOfParticles / volume;
      positions.erase(positions.begin() + selectedParticle);
    }
  }
}

void MonteCarlo::swapCBMCMove()
{
  int numberOfTrialLocations = 10;

  if (uniform() < 0.5)
  {
    insertionAttempted++;

    // trial insertion of new particle
    std::vector<double3> newParticles(numberOfTrialLocations);
    std::vector<EnergyVirial> energies(numberOfTrialLocations);
    double sumBoltzmann = 0.0;
    std::vector<double> boltzmannWeights(numberOfTrialLocations);

    // calculate energy and Boltzmann weights for all particles
    for (int i = 0; i < numberOfTrialLocations; i++)
    {
      newParticles[i] = boxSize * double3(uniform(), uniform(), uniform());
      energies[i] = particleEnergyVirial(positions, newParticles[i], numberOfParticles, boxSize, cutOff);
      boltzmannWeights[i] = std::exp(-beta * energies[i].energy);
      sumBoltzmann += boltzmannWeights[i];
    }

    if (sumBoltzmann < 1e-6) return;

    // select
    double x = uniform();
    double cumulative = 0.0;
    int selectedParticle = 0;

    for (int i = 0; i < numberOfTrialLocations; ++i)
    {
      cumulative += boltzmannWeights[i] / sumBoltzmann;
      if (x < cumulative)
      {
        selectedParticle = i;
        break;
      }
    }

    // add increase tail energy (N+1)**2 - N**2
    EnergyVirial tailEnergyDifference = (cutOffEnergyVirial(numberOfParticles + 1, boxSize, cutOffPrefactor) -
                                         cutOffEnergyVirial(numberOfParticles, boxSize, cutOffPrefactor));

    double Pacc = (pressure * beta * volume / static_cast<double>(numberOfParticles + 1.0)) * sumBoltzmann /
                  numberOfTrialLocations;
    Pacc *= std::exp(-beta * tailEnergyDifference.energy);

    if (uniform() < Pacc)
    {
      // accept
      insertionAccepted++;
      runningEnergyVirial += (tailEnergyDifference + energies[selectedParticle]);
      numberOfParticles++;
      density = numberOfParticles / volume;
      positions.push_back(newParticles[selectedParticle]);
    }
  }
  else
  {
    if (numberOfParticles == 0)
    {
      // can not delete if there are no particles
      return;
    }

    deletionAttempted++;
    // trial deletion
    int selectedParticle = static_cast<int>(uniform() * numberOfParticles);

    std::vector<double3> newParticles(numberOfTrialLocations);
    std::vector<EnergyVirial> energies(numberOfTrialLocations);
    double sumBoltzmann = 0.0;
    std::vector<double> boltzmannWeights(numberOfTrialLocations);

    newParticles[0] = positions[selectedParticle];

    // calculate energy and Boltzmann weights for all particles
    for (int i = 1; i < numberOfTrialLocations; i++)
    {
      newParticles[i] = boxSize * double3(uniform(), uniform(), uniform());
    }
    for (int i = 0; i < numberOfTrialLocations; i++)
    {
      energies[i] = particleEnergyVirial(positions, newParticles[i], numberOfParticles, boxSize, cutOff);
      boltzmannWeights[i] = std::exp(-beta * energies[i].energy);
      sumBoltzmann += boltzmannWeights[i];
    }

    // remove tail energy (N+1)**2 - N**2
    EnergyVirial tailEnergyDifference = (cutOffEnergyVirial(numberOfParticles + 1, boxSize, cutOffPrefactor) -
                                         cutOffEnergyVirial(numberOfParticles, boxSize, cutOffPrefactor));

    double Pacc = (numberOfParticles / (pressure * beta * volume)) * numberOfTrialLocations / sumBoltzmann;
    Pacc *= std::exp(-beta * tailEnergyDifference.energy);

    if (uniform() < Pacc)
    {
      deletionAccepted++;
      runningEnergyVirial += (tailEnergyDifference - energies[0]);
      numberOfParticles--;
      density = numberOfParticles / volume;
      positions.erase(positions.begin() + selectedParticle);
    }
  }
}