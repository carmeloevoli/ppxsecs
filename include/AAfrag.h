#ifndef INCLUDE_AAFRAG_H_
#define INCLUDE_AAFRAG_H_

#include "xsecs-base.h"
#include <string>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_vector.h>

class AAfrag: public XSECS {
public:
	AAfrag(Particle particle)
			: XSECS(particle) {
		if (particle == Particle::photons)
			load_table("AAfrag/spec_gamma.txt");
		else if (particle == Particle::neutrinos)
			load_table("AAfrag/spec_nutot.txt");
		else
			throw std::runtime_error("no table available for the chosen particle.");
	}

	~AAfrag() {
		gsl_vector_free(m_Ep);
		gsl_vector_free(m_Es);
		gsl_vector_free(m_sigma);
		gsl_spline2d_free(spline);
		gsl_interp_accel_free(xacc);
		gsl_interp_accel_free(yacc);
	}

	double get(double E_proj, double E_secondary) const override;

protected:
	void load_table(const std::string& filename);
	double interpolate_table(double E_proj, double E_secondary) const;

private:
	const size_t m_np = 100;
	const size_t m_ns = 1000;
	const gsl_interp2d_type *T = gsl_interp2d_bilinear;
	gsl_vector * m_Ep;
	gsl_vector * m_Es;
	gsl_vector * m_sigma;
	gsl_spline2d * spline;
	gsl_interp_accel *xacc;
	gsl_interp_accel *yacc;
};

#endif /* INCLUDE_AAFRAG_H_ */
