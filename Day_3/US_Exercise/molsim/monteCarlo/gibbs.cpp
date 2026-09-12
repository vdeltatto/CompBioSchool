#include "gibbs.h"
#include "writePDB.h"


#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <random>

void gibbsVolumeMove(MonteCarlo& systemA, MonteCarlo& systemB)
{
  systemA.volumeAttempted++;
  systemB.volumeAttempted++;

  double volumeA = systemA.volume;
  double volumeB = systemB.volume;
  double numberOfParticlesA = systemA.numberOfParticles;
  double numberOfParticlesB = systemB.numberOfParticles;
  double beta = systemA.beta;

  double totalVolume = volumeA + volumeB;
  int totalNumberOfParticles = numberOfParticlesA * numberOfParticlesB;

  double logVolumeRatio = std::log(volumeA / volumeB);
  double deltaVolume = (systemA.uniform() - 0.5) * systemA.maxVolumeChange;
  double newLogVolumeRation = logVolumeRatio + deltaVolume;
  double newVolumeRatio = std::exp(newLogVolumeRation);

  double newVolumeA = totalVolume * newVolumeRatio / (1.0 + newVolumeRatio);
  double newVolumeB = totalVolume - newVolumeA;

  if (newVolumeA <= 0.0 || newVolumeB <= 0.0) return;

  double newBoxLengthA = std::cbrt(newVolumeA);
  double newBoxLengthB = std::cbrt(newVolumeB);

  double scaleA = newBoxLengthA / systemA.boxSize;
  double scaleB = newBoxLengthB / systemB.boxSize;

  std::vector<double3> newPositionsA(systemA.positions);
  std::vector<double3> newPositionsB(systemB.positions);

  for (auto& pos : newPositionsA)
  {
    pos *= scaleA;
  }
  for (auto& pos : newPositionsB)
  {
    pos *= scaleB;
  }

  EnergyVirial newEnergyVirialA =
      systemEnergyVirial(newPositionsA, newBoxLengthA, systemA.cutOff, systemA.cutOffPrefactor);
  EnergyVirial newEnergyVirialB =
      systemEnergyVirial(newPositionsB, newBoxLengthB, systemB.cutOff, systemB.cutOffPrefactor);

  double deltaEnergyA = newEnergyVirialA.energy - systemA.runningEnergyVirial.energy;
  double deltaEnergyB = newEnergyVirialB.energy - systemB.runningEnergyVirial.energy;

  double prefactorA = (numberOfParticlesA + 1.0) * std::log(newVolumeA / volumeA);
  double prefactorB = (numberOfParticlesB + 1.0) * std::log(newVolumeB / volumeB);

  double Pacc = std::exp(prefactorA + prefactorB - beta * (deltaEnergyA + deltaEnergyB));

  if (systemA.uniform() < Pacc)
  {
    systemA.volumeAccepted++;
    systemB.volumeAccepted++;

    systemA.positions.swap(newPositionsA);
    systemB.positions.swap(newPositionsB);

    systemA.boxSize = newBoxLengthA;
    systemB.boxSize = newBoxLengthB;

    systemA.volume = newVolumeA;
    systemB.volume = newVolumeB;

    systemA.density = numberOfParticlesA / newVolumeA;
    systemB.density = numberOfParticlesB / newVolumeB;

    systemA.runningEnergyVirial = newEnergyVirialA;
    systemB.runningEnergyVirial = newEnergyVirialB;
  }
}

void gibbsSwapMove(MonteCarlo& systemA, MonteCarlo& systemB)
{
  systemA.deletionAttempted++;
  systemB.insertionAttempted++;

    double volumeA = systemA.volume;
  double volumeB = systemB.volume;
  double numberOfParticlesA = systemA.numberOfParticles;
  double numberOfParticlesB = systemB.numberOfParticles;
  double beta = systemA.beta;

  if (numberOfParticlesA == 0) return;

  int selectedParticle = static_cast<int>(systemA.uniform() * numberOfParticlesA);
  double3 oldPosition = systemA.positions[selectedParticle];
  double3 newPosition = systemB.boxSize * double3(systemA.uniform(), systemA.uniform(), systemA.uniform());

  EnergyVirial energyDifferenceA =
      particleEnergyVirial(systemA.positions, oldPosition, selectedParticle, systemA.boxSize, systemA.cutOff);
  energyDifferenceA *= -1.0;
  energyDifferenceA += (cutOffEnergyVirial(numberOfParticlesA - 1.0, systemA.boxSize, systemA.cutOffPrefactor) -
                        cutOffEnergyVirial(numberOfParticlesA, systemA.boxSize, systemA.cutOffPrefactor));

  EnergyVirial energyDifferenceB =
      particleEnergyVirial(systemB.positions, newPosition, numberOfParticlesB, systemB.boxSize, systemB.cutOff);
  energyDifferenceB += (cutOffEnergyVirial(numberOfParticlesB + 1.0, systemB.boxSize, systemB.cutOffPrefactor) -
                        cutOffEnergyVirial(numberOfParticlesB, systemB.boxSize, systemB.cutOffPrefactor));

  double deltaEnergy = energyDifferenceA.energy + energyDifferenceB.energy;
  double prefactor = numberOfParticlesA * volumeB / ((numberOfParticlesB + 1.0) * volumeB);
  double Pacc = prefactor * std::exp(-beta * deltaEnergy);

  if (systemA.uniform() < Pacc)
  {
    systemA.deletionAccepted++;
    systemB.insertionAccepted++;

    systemA.runningEnergyVirial += energyDifferenceA;
    systemB.runningEnergyVirial += energyDifferenceB;

    systemA.numberOfParticles--;
    systemB.numberOfParticles++;

    systemA.density = systemA.numberOfParticles / systemA.volume;
    systemB.density = systemB.numberOfParticles / systemB.volume;

    systemA.positions.erase(systemA.positions.begin() + selectedParticle);
    systemB.positions.push_back(newPosition);
  }
}

void runGibbsEnsemble(MonteCarlo& systemA, MonteCarlo& systemB)
{
  if (std::abs(systemA.temperature - systemB.temperature) > 1e-8)
    throw std::runtime_error("Gibbs ensemble requires both boxes at the same temperature.");

  if (std::abs(systemA.cutOff - systemB.cutOff) > 1e-8)
    throw std::runtime_error("Gibbs ensemble requires both boxes to use the same cutOff.");

  int totalNumberOfParticles = systemA.numberOfParticles + systemB.numberOfParticles;

  for (int cycle = 0; cycle < systemA.numberOfInitCycles + systemA.numberOfProdCycles; ++cycle)
  {
    systemA.cycle = cycle;
    systemB.cycle = cycle;

    for (int move = 0; move < std::max(20, totalNumberOfParticles); ++move)
    {
      double p = systemA.uniform();
      if (p < systemA.volumeProbability)
      {
        gibbsVolumeMove(systemA, systemB);
      }
      else if (p < systemA.swapProbability)
      {
        // account for move asymmetry
        if (systemA.uniform() < 0.5)
        {
          gibbsSwapMove(systemA, systemB);
        }
        else
        {
          gibbsSwapMove(systemB, systemA);
        }
      }
      else
      {
        double q = systemA.uniform() * totalNumberOfParticles;
        if (q < systemA.numberOfParticles)
        {
          systemA.translationMove();
        }
        else
        {
          systemB.translationMove();
        }
      }
    }

    if (cycle > systemA.numberOfInitCycles && cycle % 10 == 0)
    {
      systemA.computePressure();
      systemB.computePressure();
    }
    if (cycle > systemA.numberOfInitCycles)
    {
      systemA.computeChemicalPotential();
      systemB.computeChemicalPotential();
    }
    if (cycle % systemA.sampleFrequency == 0)
    {
      systemA.totalEnergyVirial =
          systemEnergyVirial(systemA.positions, systemA.boxSize, systemA.cutOff, systemA.cutOffPrefactor);
      systemB.totalEnergyVirial =
          systemEnergyVirial(systemB.positions, systemB.boxSize, systemB.cutOff, systemB.cutOffPrefactor);

      systemA.updateMoveStatistics();
      systemB.updateMoveStatistics();

      systemA.logger.debug("System A:\n" + systemA.repr());
      systemB.logger.debug("System B:\n" + systemB.repr());

      if (systemA.outputPDB)
      {
        writePDB("movie_systemA.pdb", systemA.positions, systemA.boxSize, systemA.frameNumber);
        systemA.frameNumber++;
      }

      if (systemB.outputPDB)
      {
        writePDB("movie_systemA.pdb", systemB.positions, systemB.boxSize, systemB.frameNumber);
        systemB.frameNumber++;
      }
      if (systemA.optimizeMCMoves)
      {
        systemA.optimizeMaxDisplacement();
        systemA.optimizeVolumeChange();
      }
      if (systemB.optimizeMCMoves)
      {
        systemB.optimizeMaxDisplacement();
        systemB.optimizeVolumeChange();
      }

      if (cycle > systemA.numberOfInitCycles)
      {
        systemA.driftEnergies.push_back(systemA.runningEnergyVirial.energy - systemA.totalEnergyVirial.energy);
        systemA.energies.push_back(systemA.runningEnergyVirial.energy);
        systemA.volumes.push_back(systemA.volume);
        systemA.particleCounts.push_back(systemA.numberOfParticles);
        systemA.densities.push_back(systemA.density);

        systemB.driftEnergies.push_back(systemB.runningEnergyVirial.energy - systemB.totalEnergyVirial.energy);
        systemB.energies.push_back(systemB.runningEnergyVirial.energy);
        systemB.volumes.push_back(systemB.volume);
        systemB.particleCounts.push_back(systemB.numberOfParticles);
        systemB.densities.push_back(systemB.density);
      }
    }
  }

  systemA.computePressure();
  systemB.computePressure();
  systemA.computeChemicalPotential();
  systemB.computeChemicalPotential();

  systemA.updateMoveStatistics();
  systemB.updateMoveStatistics();

  systemA.totalEnergyVirial =
      systemEnergyVirial(systemA.positions, systemA.boxSize, systemA.cutOff, systemA.cutOffPrefactor);
  systemB.totalEnergyVirial =
      systemEnergyVirial(systemB.positions, systemB.boxSize, systemB.cutOff, systemB.cutOffPrefactor);

      systemA.logger.debug("System A:\n" + systemA.repr());
      systemB.logger.debug("System B:\n" + systemB.repr());

  if (systemA.outputPDB)
  {
    writePDB("movie_systemA.pdb", systemA.positions, systemA.boxSize, systemA.frameNumber);
    systemA.frameNumber++;
  }

  if (systemB.outputPDB)
  {
    writePDB("movie_systemA.pdb", systemB.positions, systemB.boxSize, systemB.frameNumber);
    systemB.frameNumber++;
  }

  systemA.logThermodynamicalAverages();
  systemB.logThermodynamicalAverages();
}