#include "utilities.h"
#include <cmath>
#include <algorithm>

double pow2(double x) {
	return x * x;
}

/**
 * This function calculates the pp total inelastic cross section
 *
 * @param T_p is the proton kinetic energy in the LAB frame in [GeV]
 * @return inelastic cross section in [mb]
 */
double sigma_inel(double Tp) {
	double LX = std::log(Tp / Tp_th);
	double Threshold = std::max(0., 1. - std::pow(Tp_th / Tp, 1.9));
	return (Tp > Tp_th) ? (30.7 - 0.96 * LX + 0.18 * pow2(LX)) * std::pow(Threshold, 3.0) : 0;
}
