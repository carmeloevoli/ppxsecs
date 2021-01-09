// Copyright MIT license Carmelo Evoli 2020
#ifndef _KAMAE_H
#define _KAMAE_H

#include "xsecs-base.h"

class Kamae : public XSECS {
 public:
  Kamae(Particle particle) : XSECS(particle) {}

  double get(double E_proj, double E_secondary) const override;
};

#endif
