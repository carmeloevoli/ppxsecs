// Copyright MIT license Carmelo Evoli 2020
#ifndef INCLUDE_XSECS_BASE_H_
#define INCLUDE_XSECS_BASE_H_

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "Timer.h"

class XSECS {
 public:
  enum Particle { photons, neutrinos };

  XSECS(Particle particle) : m_particle(particle) {}
  virtual ~XSECS() = default;

  virtual double dsigmadE(double E_proj, double E_secondary) const = 0;

  double get(double E_proj, double E_secondary) const {
    if (E_secondary < E_proj)
      return dsigmadE(E_proj, E_secondary);
    else
      return 0;
  }

  double production_rate(double E_secondary, double alpha) const;
  void print_rate(const std::string& filename) const;
  void print(const std::string& filename) const;

 protected:
  Particle m_particle;
};

// double sigma_inel(double Tp);

#endif /* INCLUDE_XSECS_BASE_H_ */
