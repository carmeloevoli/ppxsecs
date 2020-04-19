#ifndef INCLUDE_PPGAM_H_
#define INCLUDE_PPGAM_H_

#include "xsecs-base.h"
#include <string>

class PPGam: public XSECS {
public:
	PPGam(Particle particle) :
			XSECS(particle) {
	}

	double get(double E_proj, double E_secondary) const override;

	enum InteractionModel {
		GEANT4, PYTHIA8, QGSJET, SIBYLL
	};

	void set_interaction_model(const std::string& model_name);

protected:
	InteractionModel m_intmodel = QGSJET;
};

#endif /* INCLUDE_PPGAM_H_ */
