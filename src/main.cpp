#include <fstream>
#include <iostream>
#include "AAfrag.h"
#include "kamae.h"
#include "kelner-aharonian.h"
#include "ppGam.h"

void print(const XSECS& xsecs, const std::string& filename) {
	std::ofstream outfile(filename.c_str());
	if (outfile.is_open()) {
		outfile << std::scientific << std::setprecision(4);
		outfile << "# E [GeV] - sigma [mb]\n";
		for (double E = 1; E < 3e6; E *= 1.01) {
			outfile << E << " ";
			outfile << xsecs.get(1e1, E) << " ";
			outfile << xsecs.get(1e2, E) << " ";
			outfile << xsecs.get(1e3, E) << " ";
			outfile << xsecs.get(1e5, E) << " ";
			outfile << "\n";
		}
	}
	outfile.close();
}

void printPPGam(PPGam ppgam, const std::string& filename) {
	std::ofstream outfile(filename.c_str());
	if (outfile.is_open()) {
		outfile << std::scientific << std::setprecision(4);
		outfile << "# E [GeV] - sigma [mb]\n";
		for (double E = 1; E < 3e6; E *= 1.01) {
			outfile << E << " ";
			ppgam.set_interaction_model("GEANT4");
			outfile << ppgam.get(1e1, E) << " ";
			outfile << ppgam.get(1e2, E) << " ";
			outfile << ppgam.get(1e3, E) << " ";
			outfile << ppgam.get(1e5, E) << " ";
			ppgam.set_interaction_model("PYTHIA8");
			outfile << ppgam.get(1e1, E) << " ";
			outfile << ppgam.get(1e2, E) << " ";
			outfile << ppgam.get(1e3, E) << " ";
			outfile << ppgam.get(1e5, E) << " ";
			ppgam.set_interaction_model("QGSJET");
			outfile << ppgam.get(1e1, E) << " ";
			outfile << ppgam.get(1e2, E) << " ";
			outfile << ppgam.get(1e3, E) << " ";
			outfile << ppgam.get(1e5, E) << " ";
			ppgam.set_interaction_model("SIBYLL");
			outfile << ppgam.get(1e1, E) << " ";
			outfile << ppgam.get(1e2, E) << " ";
			outfile << ppgam.get(1e3, E) << " ";
			outfile << ppgam.get(1e5, E) << " ";
			outfile << "\n";
		}
	}
	outfile.close();
}

int main() {
	try {
		Kamae kamae(Particle::photons);
		KelnerAharonian kelner(Particle::photons);
		AAfrag aafrag(Particle::photons);

		print(kamae, "specgamma-kamae.txt");
		print(kelner, "specgamma-kelner.txt");
		print(aafrag, "specgamma-aafrag.txt");

		PPGam ppgam(Particle::photons);

		printPPGam(ppgam, "specgamma-ppgam.txt");

	} catch (const std::exception& e) {
		std::cout << "\n" << e.what() << "\n";
	}
	return 0;
}
