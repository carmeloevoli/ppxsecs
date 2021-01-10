// Copyright MIT license Carmelo Evoli 2020
#ifndef INCLUDE_PPGAM_H_
#define INCLUDE_PPGAM_H_

#include <string>

#include "xsecs-base.h"

class PPGam : public XSECS {
 public:
  enum InteractionModel { GEANT4, PYTHIA8, QGSJET, SIBYLL };

  PPGam(Particle particle, InteractionModel model) : XSECS(particle), m_intmodel(model) {}

  void set_interaction_model(const std::string& model_name);

  double dsigmadE(double E_proj, double E_secondary) const override;

 protected:
  InteractionModel m_intmodel = QGSJET;
};

#endif /* INCLUDE_PPGAM_H_ */
