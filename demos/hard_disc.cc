/*
  Copyright (c) 2016-2018 Lester Hedges <lester.hedges+aabbcc@gmail.com>

  This software is provided 'as-is', without any express or implied
  warranty. In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.

  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.

  3. This notice may not be removed or altered from any source distribution.
*/

#include <chrono>
#include <fstream>
#include <iostream>

#include "MersenneTwister.h"
#include "abt/aabb_tree.hpp"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

template <class T>
using vec = abt::tree<2>::vec<T>;

/*! \file hard_disc.cpp

  An example showing the use of AABB trees for simulating the dynamics
  of a binary hard disc system where there is a large size asymmetry (10:1)
  between the particle species.
*/

// FUNCTION PROTOTYPES

// Test whether two particles overlap.
bool overlaps(abt::point<2> &,
              abt::point<2> &,
              const vec<bool> &,
              const vec<double> &,
              double);

// Compute the minimum image separation vector.
void minimumImage(vec<double> &, const vec<bool> &, const vec<double> &);

// Apply periodic boundary conditions.
void periodicBoundaries(abt::point<2> &,
                        const vec<bool> &,
                        const vec<double> &);

// Append a particle configuration to a VMD xyz file.
void printVMD(const std::string &,
              const std::vector<abt::point<2>> &,
              const std::vector<abt::point<2>> &);

// MAIN FUNCTION

int main(int argc, char **argv)
{
  using clk = std::chrono::high_resolution_clock;
  auto start = clk::now();
  // Print git commit info, if present.
#ifdef COMMIT
  std::cout << "Git commit: " << COMMIT << "\n";
#endif

  // Print git branch info, if present.
#ifdef BRANCH
  std::cout << "Git branch: " << BRANCH << "\n";
#endif

  /*****************************************************************/
  /*      Set parameters, initialise variables and objects.        */
  /*****************************************************************/

  unsigned int nSweeps = 10000;       // The number of Monte Carlo sweeps.
  unsigned int sampleInterval = 100;  // The number of sweeps per sample.
  unsigned int nSmall = 1000;         // The number of small particles.
  unsigned int nLarge = 100;          // The number of large particles.
  double diameterSmall = 1;           // The diameter of the small particles.
  double diameterLarge = 10;          // The diameter of the large particles.
  double density = 0.1;               // The system density.
  double maxDisp = 0.1;  // Maximum trial displacement (in units of diameter).

  // Total particles.
  unsigned int nParticles = nSmall + nLarge;

  // Number of samples.
  unsigned int nSamples = nSweeps / sampleInterval;

  // Particle radii.
  double radiusSmall = 0.5 * diameterSmall;
  double radiusLarge = 0.5 * diameterLarge;

  // Output formatting flag.
  unsigned int format = std::floor(std::log10(nSamples));

  // Set the periodicity of the simulation box.

  vec<bool> periodicity({true, true});

  // Work out base length of simulation box.
  double baseLength =
      std::pow((M_PI * (nSmall * diameterSmall + nLarge * diameterLarge)) /
                   (4.0 * density),
               1.0 / 2.0);
  vec<double> boxSize({baseLength, baseLength});

  // Initialise the random number generator.
  // Make it deterministic to compare across runs.
  MersenneTwister rng(0);

  // Initialise the AABB trees.
  abt::tree<2> treeSmall(maxDisp, periodicity, boxSize, nSmall);
  abt::tree<2> treeLarge(maxDisp, periodicity, boxSize, nLarge);

  // Initialise particle position vectors.
  std::vector<abt::point<2>> positionsSmall(nSmall,
                                            abt::point<2>(boxSize.size()));
  std::vector<abt::point<2>> positionsLarge(nLarge,
                                            abt::point<2>(boxSize.size()));

  /*****************************************************************/
  /*             Generate the initial AABB trees.                  */
  /*****************************************************************/

  // First the large particles.

  std::cout << "\nInserting large particles into AABB tree ...\n";
  for (unsigned int i = 0; i < nLarge; i++) {
    // Initialise the particle position vector.
    abt::point<2> position;

    // Insert the first particle directly.
    if (i == 0) {
      // Generate a random particle position.
      position[0] = boxSize[0] * rng();
      position[1] = boxSize[1] * rng();
    }

    // Check for overlaps.
    else {
      // Initialise overlap flag.
      bool isOverlap = true;

      // Keep trying until there is no overlap.
      while (isOverlap) {
        // Generate a random particle position.
        position[0] = boxSize[0] * rng();
        position[1] = boxSize[1] * rng();

        auto aabb = abt::aabb<2>::of_sphere(position, radiusLarge);

        // Query AABB overlaps.
        std::vector<unsigned int> particles = treeLarge.get_overlaps(aabb);

        // Flag as no overlap (yet).
        isOverlap = false;

        // Test overlap.
        for (unsigned int j = 0; j < particles.size(); j++) {
          // Cut-off distance.
          double cutOff = 2.0 * radiusLarge;
          cutOff *= cutOff;

          // Particles overlap.
          if (overlaps(position, positionsLarge[particles[j]], periodicity,
                       boxSize, cutOff)) {
            isOverlap = true;
            break;
          }
        }
      }
    }

    // Insert the particle into the tree.
    treeLarge.insert(i, abt::aabb2d::of_sphere(position, radiusLarge));

    // Store the position.
    positionsLarge[i] = position;
  }
  std::cout << "Tree generated!\n";

  // Now fill the gaps with small particles.

  std::cout << "\nInserting small particles into AABB tree ...\n";
  for (unsigned int i = 0; i < nSmall; i++) {
    // Initialise the particle position vector.
    abt::point<2> position;

    // Initialise overlap flag.
    bool isOverlap = true;

    // Keep trying until there is no overlap.
    while (isOverlap) {
      // Generate a random particle position.
      position[0] = boxSize[0] * rng();
      position[1] = boxSize[1] * rng();

      auto aabb = abt::aabb<2>::of_sphere(position, radiusSmall);

      // First query AABB overlaps with the large particles.
      std::vector<unsigned int> particles = treeLarge.get_overlaps(aabb);

      // Flag as no overlap (yet).
      isOverlap = false;

      // Test overlap.
      for (unsigned int j = 0; j < particles.size(); j++) {
        // Cut-off distance.
        double cutOff = radiusSmall + radiusLarge;
        cutOff *= cutOff;

        // Particles overlap.
        if (overlaps(position, positionsLarge[particles[j]], periodicity,
                     boxSize, cutOff)) {
          isOverlap = true;
          break;
        }
      }

      // Advance to next overlap test.
      if (!isOverlap) {
        // No need to test the first particle.
        if (i > 0) {
          // Now query AABB overlaps with other small particles.
          particles = treeSmall.get_overlaps(aabb);

          // Test overlap.
          for (unsigned int j = 0; j < particles.size(); j++) {
            // Cut-off distance.
            double cutOff = 2.0 * radiusSmall;
            cutOff *= cutOff;

            // Particles overlap.
            if (overlaps(position, positionsSmall[particles[j]], periodicity,
                         boxSize, cutOff)) {
              isOverlap = true;
              break;
            }
          }
        }
      }
    }

    // Insert particle into tree.
    treeSmall.insert(i, abt::aabb2d::of_sphere(position, radiusSmall));

    // Store the position.
    positionsSmall[i] = position;
  }
  std::cout << "Tree generated!\n";

  /*****************************************************************/
  /*      Perform the dynamics, updating the tree as we go.        */
  /*****************************************************************/

  // Clear the trajectory file.
  FILE *pFile;
  pFile = fopen("trajectory.xyz", "w");
  fclose(pFile);

  unsigned int sampleFlag = 0;
  unsigned int nSampled = 0;

  std::cout << "\nRunning dynamics ...\n";
  for (unsigned int i = 0; i < nSweeps; i++) {
    for (unsigned int j = 0; j < nParticles; j++) {
      // Choose a random particle.
      unsigned int particle = rng.integer(0, nParticles - 1);

      // Determine the particle type.
      unsigned int particleType = (particle < nSmall) ? 0 : 1;

      // Determine the radius of the particle.
      double radius = (particleType == 0) ? radiusSmall : radiusLarge;

      // Shift the particle index.
      if (particleType == 1)
        particle -= nSmall;

      // Initialise vectors.
      vec<double> displacement;
      abt::point<2> position;
      abt::point<2> lowerBound;
      abt::point<2> upperBound;

      // Calculate the new particle position and displacement.
      if (particleType == 0) {
        displacement[0] = maxDisp * diameterSmall * (2.0 * rng() - 1.0);
        displacement[1] = maxDisp * diameterSmall * (2.0 * rng() - 1.0);
        position[0] = positionsSmall[particle][0] + displacement[0];
        position[1] = positionsSmall[particle][1] + displacement[1];
      }
      else {
        displacement[0] = maxDisp * diameterLarge * (2.0 * rng() - 1.0);
        displacement[1] = maxDisp * diameterLarge * (2.0 * rng() - 1.0);
        position[0] = positionsLarge[particle][0] + displacement[0];
        position[1] = positionsLarge[particle][1] + displacement[1];
      }

      // Apply periodic boundary conditions.
      periodicBoundaries(position, periodicity, boxSize);

      auto aabb = abt::aabb<2>::of_sphere(position, radius);

      // Query AABB overlaps with small particles.
      std::vector<unsigned int> particles = treeSmall.get_overlaps(aabb);

      // Flag as not overlapping (yet).
      bool isOverlap = false;

      // Test overlap.
      for (unsigned int k = 0; k < particles.size(); k++) {
        // Don't test self overlap.
        if ((particleType == 1) || (particles[k] != particle)) {
          // Cut-off distance.
          double cutOff = radius + radiusSmall;
          cutOff *= cutOff;

          // Particles overlap.
          if (overlaps(position, positionsSmall[particles[k]], periodicity,
                       boxSize, cutOff)) {
            isOverlap = true;
            break;
          }
        }
      }

      // Advance to next overlap test.
      if (!isOverlap) {
        // Now query AABB overlaps with large particles.
        particles = treeLarge.get_overlaps(aabb);

        // Test overlap.
        for (unsigned int k = 0; k < particles.size(); k++) {
          // Don't test self overlap.
          if ((particleType == 0) || (particles[k] != particle)) {
            // Cut-off distance.
            double cutOff = radius + radiusLarge;
            cutOff *= cutOff;

            // Particles overlap.
            if (overlaps(position, positionsLarge[particles[k]], periodicity,
                         boxSize, cutOff)) {
              isOverlap = true;
              break;
            }
          }
        }

        // Accept the move.
        if (!isOverlap) {
          // Update the position and AABB tree.
          if (particleType == 0) {
            positionsSmall[particle] = position;
            treeSmall.update(particle, {lowerBound, upperBound});
          }
          else {
            positionsLarge[particle] = position;
            treeLarge.update(particle, {lowerBound, upperBound});
          }
        }
      }
    }

    sampleFlag++;

    if (sampleFlag == sampleInterval) {
      sampleFlag = 0;
      nSampled++;

      printVMD("trajectory.xyz", positionsSmall, positionsLarge);

      if (format == 1)
        printf("Saved configuration %2d of %2d\n", nSampled, nSamples);
      else if (format == 2)
        printf("Saved configuration %3d of %3d\n", nSampled, nSamples);
      else if (format == 3)
        printf("Saved configuration %4d of %4d\n", nSampled, nSamples);
      else if (format == 4)
        printf("Saved configuration %5d of %5d\n", nSampled, nSamples);
      else if (format == 5)
        printf("Saved configuration %6d of %6d\n", nSampled, nSamples);
    }
  }

  auto end = clk::now();
  auto elapsed_s =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count() /
      1000.;
  std::cout << "Done! Time elapsed: " << elapsed_s << "s\n";

  return (EXIT_SUCCESS);
}

// FUNCTION DEFINITIONS

bool overlaps(abt::point<2> &position1,
              abt::point<2> &position2,
              const vec<bool> &periodicity,
              const vec<double> &boxSize,
              double cutOff)
{
  // Calculate particle separation.
  vec<double> separation{(position1[0] - position2[0]),
                         (position1[1] - position2[1])};

  // Calculate minimum image separation.
  minimumImage(separation, periodicity, boxSize);

  double rSqd = separation[0] * separation[0] + separation[1] * separation[1];

  if (rSqd < cutOff)
    return true;
  else
    return false;
}

void minimumImage(vec<double> &separation,
                  const vec<bool> &periodicity,
                  const vec<double> &boxSize)
{
  for (unsigned int i = 0; i < 2; i++) {
    if (separation[i] < -0.5 * boxSize[i]) {
      separation[i] += periodicity[i] * boxSize[i];
    }
    else {
      if (separation[i] >= 0.5 * boxSize[i]) {
        separation[i] -= periodicity[i] * boxSize[i];
      }
    }
  }
}

void periodicBoundaries(abt::point<2> &position,
                        const vec<bool> &periodicity,
                        const vec<double> &boxSize)
{
  for (unsigned int i = 0; i < 2; i++) {
    if (position[i] < 0) {
      position[i] += periodicity[i] * boxSize[i];
    }
    else {
      if (position[i] >= boxSize[i]) {
        position[i] -= periodicity[i] * boxSize[i];
      }
    }
  }
}

void printVMD(const std::string &fileName,
              const std::vector<abt::point<2>> &positionsSmall,
              const std::vector<abt::point<2>> &positionsLarge)
{
  FILE *pFile;
  pFile = fopen(fileName.c_str(), "a");

  fprintf(pFile, "%lu\n\n", positionsSmall.size() + positionsLarge.size());
  for (unsigned int i = 0; i < positionsSmall.size(); i++)
    fprintf(pFile, "0 %lf %lf 0\n", positionsSmall[i][0], positionsSmall[i][1]);
  for (unsigned int i = 0; i < positionsLarge.size(); i++)
    fprintf(pFile, "1 %lf %lf 0\n", positionsLarge[i][0], positionsLarge[i][1]);

  fclose(pFile);
}
