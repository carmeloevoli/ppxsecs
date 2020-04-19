#ifndef INCLUDE_UTILITIES_H_
#define INCLUDE_UTILITIES_H_

// definition of some constants

static const double m_p = 0.938272; // GeV, proton mass (taken from PDG)
static const double m_pi = 0.134976; // GeV, pi0 mass (taken from PDG)
static const double Tp_th = 2.0 * m_pi + (m_pi * m_pi) / (2.0 * m_p); // in GeV, the proton threshold energy in the LAB frame

// useful functions

double pow2(double x);
double sigma_inel(double Tp);

#endif /* INCLUDE_UTILITIES_H_ */
