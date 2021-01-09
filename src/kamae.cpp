// Copyright Carmelo Evoli 2020
#include "kamae.h"

#include "cparamlib.h"

/**
 * Calculates gamma-ray differential cross
 * section as a function of the proton kinetic energy.
 *
 * @param E_proj  Proton kinetic energy in GeV.
 * @param E_gamma Gamma-ray energy in GeV.
 * @return Cross section dsigma/dE in mb/GeV
 */
double gamma_xsec(double E_proj, double E_gamma) {
  PARTICLE_IDS par = ID_GAMMA;
  PARAMSET parameters;
  gamma_param(E_proj, &parameters);
  return sigma_incl_tot(par, E_gamma, E_proj, &parameters) / E_gamma;
}

/**
 * Calculates total neutrinos differential cross
 * section as a function of the proton kinetic energy.
 *
 * @param E_proj Proton kinetic energy in GeV.
 * @param E_nu   Neutrino energy in GeV.
 * @return Cross section dsigma/dE in mb/GeV
 */
double neutrinos_xsec(double E_proj, double E_nu) {
  double value = 0;
  PARAMSET parameters;
  nue_param(E_proj, &parameters);
  value += sigma_incl_tot(ID_NUE, E_nu, E_proj, &parameters);
  numu_param(E_proj, &parameters);
  value += sigma_incl_tot(ID_NUMU, E_nu, E_proj, &parameters);
  antinue_param(E_proj, &parameters);
  value += sigma_incl_tot(ID_ANTINUE, E_nu, E_proj, &parameters);
  antinumu_param(E_proj, &parameters);
  value += sigma_incl_tot(ID_ANTINUMU, E_nu, E_proj, &parameters);
  return value / E_nu;
}

double Kamae::get(double E_proj, double E_secondary) const {
  if (E_secondary > E_proj) return 0;

  double value = 0;

  if (m_particle == Particle::photons) {
    value = gamma_xsec(E_proj, E_secondary);
  } else if (m_particle == Particle::neutrinos) {
    value = neutrinos_xsec(E_proj, E_secondary);
  }
  return value;
}
