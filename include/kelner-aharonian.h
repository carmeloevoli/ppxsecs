#ifndef INCLUDE_KELNER_AHARONIAN_H_
#define INCLUDE_KELNER_AHARONIAN_H_

#include "xsecs-base.h"

class KelnerAharonian: public XSECS {
public:
	KelnerAharonian(Particle particle)
			: XSECS(particle) {
	}

	double get(double E_proj, double E_secondary) const override;
};

#endif /* INCLUDE_KELNER_AHARONIAN_H_ */
