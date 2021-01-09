// Copyright MIT license Carmelo Evoli 2020
#ifndef INCLUDE_XSECS_BASE_H_
#define INCLUDE_XSECS_BASE_H_

enum Particle { photons, neutrinos };

class XSECS {
 protected:
  Particle m_particle;

 public:
  XSECS(Particle particle) : m_particle(particle) {}

  virtual ~XSECS() = default;

  virtual double get(double E_proj, double E_secondary) const = 0;
};

double sigma_inel(double Tp);

#endif /* INCLUDE_XSECS_BASE_H_ */
