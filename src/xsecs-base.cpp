#include "xsecs-base.h"

double XSECS::production_rate(double E_secondary, double alpha) const {
  double value = 0;
  for (double E_p = E_secondary; E_p < 1e6; E_p *= 1.01) {
    const auto sigma = get(E_p, E_secondary);
    if (sigma != sigma) std::cout << sigma << " " << E_p << " " << E_secondary << "\n";
    value += std::pow(E_p, 1. - alpha) * sigma;
  }
  return value;
}

void XSECS::print_rate(const std::string& filename) const {
  std::ofstream outfile(filename.c_str());
  if (outfile.is_open()) {
    outfile << std::scientific << std::setprecision(4);
    outfile << "# alpha - rate [au]\n";
    for (double alpha = 2.0; alpha <= 3.0; alpha += 0.01) {
      outfile << alpha << " ";
      outfile << production_rate(1e3, alpha) << " ";
      outfile << production_rate(1e4, alpha) << " ";
      outfile << production_rate(1e5, alpha) << " ";
      outfile << "\n";
    }
  }
  outfile.close();
}

void XSECS::print(const std::string& filename) const {
  Timer timer;
  std::ofstream outfile(filename.c_str());
  if (outfile.is_open()) {
    outfile << std::scientific << std::setprecision(4);
    outfile << "# E [GeV] - sigma [mb/GeV]\n";
    for (double E = 0.01; E < 3e6; E *= 1.01) {
      outfile << E << " ";
      outfile << get(1e1, E) << " ";
      outfile << get(1e2, E) << " ";
      outfile << get(1e3, E) << " ";
      outfile << get(1e4, E) << " ";
      outfile << get(1e5, E) << " ";
      outfile << "\n";
    }
  }
  outfile.close();
}
