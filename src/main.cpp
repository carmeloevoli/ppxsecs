// Copyright MIT license Carmelo Evoli 2020
#include "AAfrag.h"
#include "kamae.h"
#include "kelner-aharonian.h"
#include "ppGam.h"

int main() {
  try {
    {
      Kamae kamae(XSECS::Particle::photons);
      kamae.print("specgamma-kamae.txt");
      kamae.print_rate("rategamma-kamae.txt");
    }
    {
      Kamae kamae(XSECS::Particle::neutrinos);
      kamae.print("specnu-kamae.txt");
      kamae.print_rate("ratenu-kamae.txt");
    }
    {
      KelnerAharonian kelner(XSECS::Particle::photons);
      kelner.print("specgamma-kelner.txt");
      kelner.print_rate("rategamma-kelner.txt");
    }
    {
      KelnerAharonian kelner(XSECS::Particle::neutrinos);
      kelner.print("specnu-kelner.txt");
      kelner.print_rate("ratenu-kelner.txt");
    }
    {
      AAfrag aafrag(XSECS::Particle::photons);
      aafrag.print("specgamma-aafrag.txt");
      aafrag.print_rate("rategamma-aafrag.txt");
    }
    {
      AAfrag aafrag(XSECS::Particle::neutrinos);
      aafrag.print("specnu-aafrag.txt");
      aafrag.print_rate("ratenu-aafrag.txt");
    }
    {
      PPGam ppgam(XSECS::Particle::photons, PPGam::InteractionModel::GEANT4);
      ppgam.print("specgamma-ppgam-GEANT4.txt");
    }
    {
      PPGam ppgam(XSECS::Particle::neutrinos, PPGam::InteractionModel::GEANT4);
      ppgam.print("specnu-ppgam-GEANT4.txt");
    }
    {
      PPGam ppgam(XSECS::Particle::photons, PPGam::InteractionModel::PYTHIA8);
      ppgam.print("specgamma-ppgam-PYTHIA8.txt");
    }
    {
      PPGam ppgam(XSECS::Particle::neutrinos, PPGam::InteractionModel::PYTHIA8);
      ppgam.print("specnu-ppgam-PYTHIA8.txt");
    }
    {
      PPGam ppgam(XSECS::Particle::photons, PPGam::InteractionModel::QGSJET);
      ppgam.print("specgamma-ppgam-QGSJET.txt");
    }
    {
      PPGam ppgam(XSECS::Particle::neutrinos, PPGam::InteractionModel::QGSJET);
      ppgam.print("specnu-ppgam-QGSJET.txt");
    }
    {
      PPGam ppgam(XSECS::Particle::photons, PPGam::InteractionModel::SIBYLL);
      ppgam.print("specgamma-ppgam-SIBYLL.txt");
    }
    {
      PPGam ppgam(XSECS::Particle::neutrinos, PPGam::InteractionModel::SIBYLL);
      ppgam.print("specnu-ppgam-SIBYLL.txt");
    }
  } catch (const std::exception& e) {
    std::cout << "\n" << e.what() << "\n";
  }
  return 0;
}
