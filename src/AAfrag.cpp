// Copyright Carmelo Evoli 2020
#include "AAfrag.h"

#include <cmath>
#include <fstream>
#include <iostream>

bool file_exists(const std::string& name) {
  std::ifstream f(name.c_str());
  return f.good();
}

void AAfrag::load_table(const std::string& filename) {
  m_Ep = gsl_vector_alloc(m_np);
  m_Es = gsl_vector_alloc(m_ns);

  const double E_min = 1e-3;

  for (size_t i = 0; i < m_np; i++) gsl_vector_set(m_Ep, i, 4e3 * E_min * std::pow(1.2, i));

  for (size_t j = 0; j < m_ns; j++) gsl_vector_set(m_Es, j, E_min * std::pow(1.025, j));

  m_sigma = gsl_vector_alloc(m_np * m_ns);
  spline = gsl_spline2d_alloc(T, m_np, m_ns);
  xacc = gsl_interp_accel_alloc();
  yacc = gsl_interp_accel_alloc();

  if (file_exists(filename)) {
    std::fstream file_to_read(filename);
    double E_p, E_s, sigma;
    for (size_t i = 0; i < m_np; i++) {
      for (size_t j = 0; j < m_ns; j++) {
        file_to_read >> E_p >> E_s >> sigma;
        gsl_spline2d_set(spline, m_sigma->data, i, j, sigma);
      }
    }
    gsl_spline2d_init(spline, m_Ep->data, m_Es->data, m_sigma->data, m_np, m_ns);
  } else {
    throw std::runtime_error("file for reading AAfrag model cannot be found.");
  }
}

double AAfrag::interpolate_table(double E_proj, double E_secondary) const {
  if (E_secondary < 1e-3 || E_secondary > 0.51658e8)
    throw std::runtime_error("E_secondary outside implemented energy range");
  if (E_proj < 0.40000e+01 || E_proj > 0.27606e+09) throw std::runtime_error("E_proj outside implemented energy range");
  return gsl_spline2d_eval(spline, E_proj, E_secondary, xacc, yacc);
}

double AAfrag::get(double E_proj, double E_secondary) const {
  if (E_secondary > E_proj) return 0;
  return interpolate_table(E_proj, E_secondary) / E_secondary;
}
